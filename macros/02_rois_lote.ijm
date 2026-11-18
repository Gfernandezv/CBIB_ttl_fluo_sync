// ============================================================================
// 02_rois_lote.ijm  -  RoiSet y ROI_Results de TODOS los stacks limpios de una carpeta
//                      y sus subcarpetas (ImageJ/Fiji)
//
// Para cada archivo *_clean_quant.tif (salida de 01_limpieza_lote.ijm o del notebook de Python):
// Esta macro RETOMA lo que dejo el paso 1 (01_limpieza_lote.ijm): por defecto procesa los
// *_clean_quant.tif que estan en las carpetas "limpieza_ij" (INPUT_SUBDIR) y usa la vista
// <nombre>_clean_view.tif que dejo junto a cada uno. No necesita nada mas.
//   1. Obtiene una imagen de etiquetas (un numero por celula), segun ROI_SOURCE:
//        "auto"      : (por defecto) usa <nombre>_rois_labels.tif si existe junto al stack; si no,
//                      segmenta con StarDist la vista que dejo el paso 1.
//        "etiquetas" : abre <nombre>_rois_labels.tif, que debe estar junto al stack
//                      (p. ej. el que deja el notebook de Python con StarDist). VALIDADO.
//        "stardist"  : segmenta con el plugin StarDist de Fiji la media de <nombre>_clean_view.tif.
//                      NO VALIDADO: puede fallar segun la instalacion y, si lanza una excepcion,
//                      detiene el lote (se puede reanudar).
//        "existente" : no crea ROI; usa un RoiSet que ya existe junto al stack (RoiSet.zip o
//                      <nombre>_RoiSet.zip) y solo mide.
//   2. Convierte las etiquetas en ROI: descarta objetos de area menor que MIN_AREA, deja la
//      componente mayor de cada uno, traza el contorno por los bordes de pixel (el ROI reproduce
//      la mascara exactamente) y los ordena de arriba hacia abajo por el centroide. La posicion k
//      en el ROI Manager es ROIk.
//   3. Guarda el RoiSet y, si DO_MEASURE, mide Area, IntDen y RawIntDen por ROI y frame con el
//      mismo formato que el Multi Measure de ImageJ (ROI_Results.csv).
//
// Salidas, en <muestra>/<OUT_SUBDIR>/ (por defecto rois_ij), al lado de <muestra>/limpieza_ij/:
//   <nombre>_RoiSet.zip           (o RoiSet.zip si FIXED_NAMES = true)
//   <nombre>_ROI_Results.csv      (o ROI_Results.csv si FIXED_NAMES = true)
//   <nombre>_rois_table_ij.csv    centroide y area de cada ROI
// y en la carpeta raiz: rois_lote_log.csv, un DIARIO con fecha (una fila por stack evaluado).
//
// IMPORTANTE: los ROI estan en las coordenadas del stack alineado y recortado (*_clean_quant.tif).
// Ademas, StarDist define candidatos por forma e intensidad, no por respuesta: la inclusion
// manual sigue siendo tuya.
//
// CONECTIVIDAD: para quedarse con la componente mayor de cada objeto, Analyze Particles usa conectividad
// 8 (un pixel unido solo en diagonal cuenta como parte del objeto). El modulo de Python (stardist_rois.py)
// usa conectividad 4 y lo descarta. En un stack real, 1 de 278 ROI difirio en 1 pixel por este motivo
// (el resto fue identico); si ocurre, ese ROI tiene 1 pixel mas en ImageJ y su fila de ROI_Results difiere.
//
// STARDIST: la macro comprueba al empezar si el comando 'StarDist 2D' esta instalado. Si el plugin
// lanza una excepcion, la macro se detiene (no se puede capturar); antes de llamarlo deja
// "stardist_estado.txt" con INICIADO, y la ejecucion siguiente lo detecta, no vuelve a llamar al plugin
// y sigue con los stacks que tengan etiquetas. Borra ese archivo cuando hayas resuelto el problema.
//
// ROBUSTEZ: una macro no puede capturar errores. Un archivo vacio se omite con un aviso. Un archivo
// corrupto, o un error grave (p. ej. falta de memoria), detiene la macro; pero antes de procesar cada
// stack se deja en el log una fila provisional INTERRUMPIDO. Al volver a ejecutar, ese stack se salta
// (y se explica por que) y se continua con el resto, en vez de bloquearse siempre en el mismo.
// Para reintentarlo, borra su fila INTERRUMPIDO del log.
//
// QUE SE SALTA (y por que):
//   - Carpetas de EXCLUDE_DIRS (por nombre, sin distinguir mayusculas).
//   - Stacks que ya tienen resultados, salvo OVERWRITE = true (se puede interrumpir y reanudar).
//   - Stacks sin etiquetas / sin vista / sin RoiSet segun el modo, o con etiquetas de otro tamano.
//
// Uso interactivo: ejecutar y elegir la carpeta raiz. Si el selector da problemas, escribe la
//                  ruta en ROOT_PATH. Si no encuentra stacks, el Log (Window > Log) explica que vio.
// Uso por lotes  : argumento "raiz|fuente|medir" (los dos ultimos opcionales; medir = 1 | 0).
//
// DISENO (importante): ninguna funcion de esta macro escribe variables globales. Los resultados se
// devuelven como texto y los arreglos se pasan como parametro. Motivo: segun el camino por el que se
// ejecute la macro, las declaraciones "var" no se registran como globales; entonces una asignacion
// dentro de una funcion, como  lista = Array.concat(...),  crea una variable LOCAL y la global queda
// vacia, mientras que los incrementos (n++) si funcionan. (Esto hizo fallar en silencio a una
// version anterior de 01_limpieza_lote.ijm.)
// ============================================================================

// ---------- Parametros ----------
var ROI_SOURCE       = "auto";        // "auto" | "etiquetas" | "stardist" | "existente"
                                      //   auto = usa <nombre>_rois_labels.tif si existe; si no, segmenta con StarDist
var INPUT_SUBDIR     = "limpieza_ij"; // carpeta donde el paso 1 (01_limpieza_lote.ijm) deja sus stacks; solo se
                                      // procesan los *_clean_quant.tif que estan en ella. "" = en cualquier carpeta
var MIN_AREA         = 0;             // area minima de un ROI en px; 0 = automatico: max(20, 25 % de la mediana)
var DO_MEASURE       = true;          // medir y guardar ROI_Results.csv
var OUT_SUBDIR       = "rois_ij";
var FIXED_NAMES      = false;         // true = RoiSet.zip y ROI_Results.csv (solo si hay un stack por carpeta)
var EXCLUDE_DIRS     = "ROIs,batch_analysis,rois_ij";   // nombres separados por coma
var OVERWRITE        = false;
var ROOT_PATH        = "";            // ruta de la carpeta raiz; vacia = abrir el selector. Se puede pegar con barras invertidas
var VERBOSE          = true;          // imprime cada carpeta visitada en el Log
// Solo para ROI_SOURCE = "stardist" (valores del modelo 2D_versatile_fluo; coinciden con Python)
var PROJECTION       = "mean";        // resumen temporal de la vista: "mean" | "max"
var STARDIST_MODEL   = "Versatile (fluorescent nuclei)";
var PROB_THRESH      = 0.479071;
var NMS_THRESH       = 0.3;
var EXCLUDE_BOUNDARY = 0;             // 0 = no excluir objetos que tocan el borde (como en Python)

// ---------- Robustez ante archivos que detienen la macro ----------
// Una macro no puede capturar el error de open() ante un archivo corrupto: se detiene. Para que eso no
// bloquee las ejecuciones siguientes, antes de procesar cada archivo se deja en el log una fila
// provisional "INTERRUMPIDO"; si la macro termina bien con el archivo, el log se reescribe sin ella.
function ahora() {
    getDateAndTime(yr, mo, dw, dd, hh, mi, ss, ms);
    return "" + yr + "-" + IJ.pad(mo + 1, 2) + "-" + IJ.pad(dd, 2) + " " + IJ.pad(hh, 2) + ":" + IJ.pad(mi, 2);
}

// true si la ultima fila del log para este archivo es un marcador INTERRUMPIDO (la macro se detuvo con el)
function ultimaFilaInterrumpida(logTxt, rel) {
    marca = "\"" + rel + "\"";
    i = lastIndexOf(logTxt, marca);
    if (i < 0) return false;
    resto = substring(logTxt, i + lengthOf(marca), minOf(lengthOf(logTxt), i + lengthOf(marca) + 16));
    return startsWith(resto, ",\"INTERRUMPIDO");
}

// ---------- Funciones auxiliares ----------
// mediana de las entradas positivas cnt[desde..hasta]
function medianaPositivos(cnt, desde, hasta) {
    pos = newArray(0);
    for (i = desde; i <= hasta; i++) if (cnt[i] > 0) pos = Array.concat(pos, cnt[i]);
    if (pos.length == 0) return 0;
    Array.sort(pos);
    n = pos.length;
    if (n % 2 == 1) return pos[(n - 1) / 2];
    return (pos[n / 2 - 1] + pos[n / 2]) / 2;
}

function carpetaExcluida(nombre) {
    if (endsWith(nombre, "/") || endsWith(nombre, "\\")) nombre = substring(nombre, 0, lengthOf(nombre) - 1);
    ex = split(EXCLUDE_DIRS, ",");
    for (i = 0; i < ex.length; i++) if (toLowerCase(String.trim(ex[i])) == toLowerCase(nombre)) return true;
    return false;
}

// true si dir es la carpeta de entrada (la que dejo el paso 1), o si no se restringe
function carpetaEntrada(dir) {
    if (INPUT_SUBDIR == "") return true;
    partes = split(dir, "/");
    if (partes.length == 0) return false;
    // (no escribir "return texto == texto": esta comparacion devuelve el primer texto en lugar de true/false)
    if (toLowerCase(partes[partes.length - 1]) == toLowerCase(INPUT_SUBDIR)) return true;
    return false;
}

function extensionDe(nm) {
    pt = lastIndexOf(nm, ".");
    if (pt < 0) return "(sin extension)";
    return toLowerCase(substring(nm, pt));
}

// Devuelve TODO lo que vio, como texto (una linea por hecho, campos separados por tabulador):
//   Q <ruta>                      stack *_clean_quant.tif a procesar
//   N <ruta>                      stack *_clean_quant.tif ignorado por no estar en INPUT_SUBDIR
//   X <ruta>                      carpeta excluida
//   E <extension>                 un archivo (para contar extensiones)
//   D <ruta> <arch> <sub> <quant> carpeta visitada y sus conteos
// Solo usa variables locales y valores devueltos: no depende de variables globales.
function buscarStacks(dir) {
    nombres = getFileList(dir);
    acum = "";
    nArch = 0; nSub = 0; nQ = 0;
    for (i = 0; i < nombres.length; i++) {
        nm = nombres[i];
        if (endsWith(nm, "/") || endsWith(nm, "\\")) {
            if (carpetaExcluida(nm)) {
                acum = acum + "X\t" + dir + nm + "\n";
            } else {
                nSub++;
                acum = acum + buscarStacks(dir + nm);
            }
        } else {
            nArch++;
            acum = acum + "E\t" + extensionDe(nm) + "\n";
            bajo = toLowerCase(nm);
            if (endsWith(bajo, "_clean_quant.tif") || endsWith(bajo, "_clean_quant.tiff")) {
                if (carpetaEntrada(dir)) {
                    acum = acum + "Q\t" + dir + nm + "\n";
                    nQ++;
                } else {
                    acum = acum + "N\t" + dir + nm + "\n";        // es un stack limpio, pero no esta en INPUT_SUBDIR
                }
            }
        }
    }
    acum = acum + "D\t" + dir + "\t" + nArch + "\t" + nSub + "\t" + nQ + "\n";
    return acum;
}

// ---------- Etiquetas -> ROI ----------
// Requiere una imagen abierta con el titulo "labels" (2D, 8 o 16 bits). Deja los ROI en el ROI Manager,
// guarda el RoiSet y una tabla de centroide y area. Devuelve "OK|objetos|rois|descartados|area_minima"
// o el motivo por el que se omitio.
function etiquetasARois(outZip, outTabla) {
    selectImage("labels");
    getDimensions(lw, lh, lc, ls, lf);
    if (ls * lf * lc > 1) return "omitido: la imagen de etiquetas debe ser 2D";
    bd = bitDepth();
    if (bd != 8 && bd != 16) return "omitido: etiquetas de " + bd + " bits (se esperan 8 o 16)";
    getStatistics(area, mean, mn, nEtiquetas);
    if (nEtiquetas < 1) return "omitido: la imagen de etiquetas no tiene objetos";
    if (bd == 16) getHistogram(vals, cnt, 65536);
    else getHistogram(vals, cnt, 256);                // cnt[k] = area (px) de la etiqueta k

    if (MIN_AREA > 0) minArea = MIN_AREA;
    else minArea = maxOf(20, 0.25 * medianaPositivos(cnt, 1, nEtiquetas));

    roiManager("reset");
    run("Set Measurements...", "area centroid redirect=None decimal=3");
    descartados = 0;
    for (k = 1; k <= nEtiquetas; k++) {
        selectImage("labels");
        run("Select None");                           // sin esto, analizar se restringe al ROI anterior
        areaK = cnt[k];
        if (areaK == 0) continue;                     // etiqueta vacia
        if (areaK < minArea) { descartados++; continue; }

        setThreshold(k, k);
        n0 = roiManager("count");
        run("Analyze Particles...", "size=1-Infinity pixel include add");   // un ROI por componente conexa
        resetThreshold();
        n1 = roiManager("count");
        if (n1 == n0) continue;

        // quedarse con la componente de mayor area
        mejor = n0; mejorA = -1;
        for (j = n0; j < n1; j++) {
            roiManager("select", j);
            getStatistics(a);
            if (a > mejorA) { mejorA = a; mejor = j; }
        }
        for (j = n1 - 1; j > mejor; j--) { roiManager("select", j); roiManager("delete"); }
        for (j = mejor - 1; j >= n0; j--) { roiManager("select", j); roiManager("delete"); }
        idx = roiManager("count") - 1;

        // clave de orden: centroide y, luego x (rellenada con ceros para ordenar como texto)
        roiManager("select", idx);
        List.setMeasurements;
        cx = List.getValue("X"); cy = List.getValue("Y");
        roiManager("rename", IJ.pad(round(cy * 1000), 9) + "_" + IJ.pad(round(cx * 1000), 9));
    }
    nRoi = roiManager("count");
    if (nRoi == 0) return "omitido: no quedo ningun ROI (area minima " + d2s(minArea, 1) + " px)";

    roiManager("sort");
    for (j = 0; j < nRoi; j++) {                      // nombre estilo ImageJ "yyyy-xxxx" (centro del cuadro)
        roiManager("select", j);
        getSelectionBounds(bx, by, bw, bh);
        roiManager("rename", IJ.pad(floor(by + bh / 2), 4) + "-" + IJ.pad(floor(bx + bw / 2), 4));
    }
    roiManager("deselect");
    roiManager("save", outZip);

    // tabla de centroide y area, ROI por ROI (roiManager("measure") desde una macro solo devolvia el primero)
    selectImage("labels");
    run("Clear Results");
    for (j = 0; j < nRoi; j++) {
        roiManager("select", j);
        List.setMeasurements;
        setResult("Area", j, List.getValue("Area"));
        setResult("X", j, List.getValue("X"));
        setResult("Y", j, List.getValue("Y"));
    }
    run("Select None");
    updateResults();
    saveAs("Results", outTabla);
    return "OK|" + nEtiquetas + "|" + nRoi + "|" + descartados + "|" + d2s(minArea, 1);
}

// ---------- Medicion ----------
// Requiere el stack cuantitativo abierto con el titulo "quant". Devuelve "OK|frames|rois" o el motivo.
function medirRois(zipPath, outCsv) {
    selectImage("quant");
    n = nSlices;
    roiManager("reset");
    roiManager("open", zipPath);
    nRoi = roiManager("count");
    if (nRoi == 0) return "omitido: el RoiSet no tiene ROI";
    suma = newArray(nRoi * n);            // suma[k * n + (frame - 1)] = suma de valores del ROI k (0-based)
    area = newArray(nRoi);
    selectImage("quant");
    for (k = 0; k < nRoi; k++) {
        roiManager("select", k);          // el ROI se mantiene al cambiar de slice
        for (i = 1; i <= n; i++) {
            setSlice(i);
            getRawStatistics(nPix, media);
            if (i == 1) area[k] = nPix;
            suma[k * n + (i - 1)] = round(nPix * media);    // suma exacta de enteros
        }
    }
    run("Select None");
    run("Clear Results");
    for (i = 0; i < n; i++) {
        for (k = 0; k < nRoi; k++) {
            v = suma[k * n + i];
            setResult("Area" + (k + 1), i, area[k]);
            setResult("IntDen" + (k + 1), i, v);
            setResult("RawIntDen" + (k + 1), i, v);
        }
    }
    updateResults();
    saveAs("Results", outCsv);
    return "OK|" + n + "|" + nRoi;
}

// ---------- Procesar un stack ----------
// Devuelve "OK|objetos|rois|descartados|area_minima|frames_medidos|segundos" o el motivo de omision.
// No usa exit(): un stack problematico no detiene el lote.
function procesarUno(qPath, hayStarDist, testigo) {
    t0 = getTime();
    dir = File.getDirectory(qPath);
    nombre = File.getName(qPath);
    stem = substring(nombre, 0, indexOf(toLowerCase(nombre), "_clean_quant"));
    // Las salidas van al lado de la carpeta del paso 1 (<muestra>/limpieza_ij y <muestra>/rois_ij), no dentro de ella
    base = dir;
    if (INPUT_SUBDIR != "") {
        t = dir;
        if (endsWith(t, "/") || endsWith(t, "\\")) t = substring(t, 0, lengthOf(t) - 1);
        i1 = lastIndexOf(t, "/");
        i2 = lastIndexOf(t, "\\");
        if (i2 > i1) i1 = i2;
        if (i1 >= 0) {
            if (toLowerCase(substring(t, i1 + 1)) == toLowerCase(INPUT_SUBDIR)) base = substring(t, 0, i1 + 1);
        }
    }
    outDir = base + OUT_SUBDIR + "/";
    if (FIXED_NAMES) { zipOut = outDir + "RoiSet.zip"; csvOut = outDir + "ROI_Results.csv"; }
    else { zipOut = outDir + stem + "_RoiSet.zip"; csvOut = outDir + stem + "_ROI_Results.csv"; }
    tablaOut = outDir + stem + "_rois_table_ij.csv";

    if (!OVERWRITE) {
        if (ROI_SOURCE == "existente") {
            if (File.exists(csvOut)) return "omitido: ya existe " + OUT_SUBDIR + "/" + File.getName(csvOut);
        } else if (File.exists(zipOut) && (!DO_MEASURE || File.exists(csvOut))) {
            return "omitido: ya existe " + OUT_SUBDIR + "/" + File.getName(zipOut);
        }
    }

    // comprobaciones previas, antes de abrir nada pesado
    if (File.length(qPath) == 0) return "omitido: el stack esta vacio (0 bytes)";
    labelsPath = dir + stem + "_rois_labels.tif";
    viewPath = dir + stem + "_clean_view.tif";
    zipFinal = zipOut;
    fuente = ROI_SOURCE;
    if (fuente == "auto") {
        if (File.exists(labelsPath)) fuente = "etiquetas";
        else fuente = "stardist";
    }
    if (fuente == "etiquetas") {
        if (!File.exists(labelsPath)) return "omitido: falta " + stem + "_rois_labels.tif (etiquetas de StarDist, p. ej. del notebook de Python)";
        if (File.length(labelsPath) == 0) return "omitido: " + stem + "_rois_labels.tif esta vacio (0 bytes)";
    } else if (fuente == "stardist") {
        if (!File.exists(viewPath)) return "omitido: no hay " + stem + "_rois_labels.tif ni " + stem + "_clean_view.tif (se necesita la vista para segmentar con StarDist)";
        if (File.length(viewPath) == 0) return "omitido: " + stem + "_clean_view.tif esta vacio (0 bytes)";
        if (!hayStarDist) return "omitido: no hay " + stem + "_rois_labels.tif y StarDist no esta disponible (no instalado o fallo previo; ver el inicio del Log)";
    } else if (fuente == "existente") {
        if (File.exists(dir + "RoiSet.zip")) zipFinal = dir + "RoiSet.zip";
        else if (File.exists(dir + stem + "_RoiSet.zip")) zipFinal = dir + stem + "_RoiSet.zip";
        else if (File.exists(zipOut)) zipFinal = zipOut;
        else if (File.exists(base + "RoiSet.zip")) zipFinal = base + "RoiSet.zip";
        else return "omitido: no hay RoiSet.zip ni " + stem + "_RoiSet.zip junto al stack";
    } else {
        return "omitido: ROI_SOURCE desconocido (" + ROI_SOURCE + ")";
    }
    File.makeDirectory(outDir);

    open(qPath);
    rename("quant");
    if (bitDepth() == 32) return "omitido: el stack cuantitativo es de 32 bits (se esperan enteros)";
    qw = getWidth(); qh = getHeight();
    if (nSlices < 2) return "omitido: no es un stack";

    nDet = 0; nRoi = 0; descartados = 0; minAreaUsada = "0";
    if (fuente != "existente") {
        if (fuente == "etiquetas") {
            open(labelsPath);
            rename("labels");
        } else {
            // plugin StarDist de Fiji sobre la media de la vista que dejo el paso 1 (la llamada al plugin NO esta validada)
            open(viewPath);
            rename("vista");
            nv = nSlices;
            if (PROJECTION == "max") {
                run("Z Project...", "projection=[Max Intensity]");
            } else {
                run("Z Project...", "projection=[Sum Slices]");
                run("Divide...", "value=" + nv);
            }
            rename("media2d");
            close("vista");
            File.saveString("INICIADO " + ahora() + "\n", testigo);     // si la macro se detiene en la llamada, queda constancia
            run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], "
                + "args=['input':'media2d', 'modelChoice':'" + STARDIST_MODEL + "', "
                + "'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', "
                + "'probThresh':'" + PROB_THRESH + "', 'nmsThresh':'" + NMS_THRESH + "', "
                + "'outputType':'Label Image', 'nTiles':'1', 'excludeBoundary':'" + EXCLUDE_BOUNDARY + "', "
                + "'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', "
                + "'showProbAndDist':'false'], process=[false]");
            File.saveString("OK " + ahora() + "\n", testigo);           // el plugin volvio: no detuvo la macro
            selectImage("Label Image");
            rename("labels");
        }
        selectImage("labels");
        if (getWidth() != qw || getHeight() != qh)
            return "omitido: las etiquetas miden " + getWidth() + "x" + getHeight() + " y el stack " + qw + "x" + qh + " (no corresponden)";
        r = etiquetasARois(zipOut, tablaOut);
        if (!startsWith(r, "OK|")) return r;
        c = split(r, "|");
        nDet = c[1]; nRoi = c[2]; descartados = c[3]; minAreaUsada = c[4];
        close("labels");
    }

    nFramesMed = 0;
    if (DO_MEASURE) {
        r2 = medirRois(zipFinal, csvOut);
        if (!startsWith(r2, "OK|")) return r2;
        c2 = split(r2, "|");
        nFramesMed = c2[1];
        if (nRoi == 0) nRoi = c2[2];
    }
    seg = round((getTime() - t0) / 1000);
    return "OK|" + nDet + "|" + nRoi + "|" + descartados + "|" + minAreaUsada + "|" + nFramesMed + "|" + seg + "|" + fuente;
}

// ---------- Programa principal ----------
arg = getArgument();
if (arg != "") {
    p = split(arg, "|");
    raiz = p[0];
    if (p.length > 1) ROI_SOURCE = p[1];
    if (p.length > 2) DO_MEASURE = (p[2] == "1");
} else if (ROOT_PATH != "") {
    raiz = ROOT_PATH;
} else {
    raiz = getDirectory("Elige la carpeta raiz (se buscan stacks *_clean_quant.tif en ella y en sus subcarpetas)");
}
print("Carpeta recibida: [" + raiz + "]");
if (raiz == "" || raiz == "null") exit("No se eligio ninguna carpeta (el selector devolvio un valor vacio).\nDefine ROOT_PATH al inicio de la macro con la ruta de la carpeta.");
raiz = replace(raiz, "\\", "/");                 // ImageJ acepta "/" tambien en Windows; asi las rutas son uniformes
if (!endsWith(raiz, "/")) raiz = raiz + "/";
if (!File.isDirectory(raiz)) exit("No existe la carpeta: " + raiz + "\nRevisa la ruta (o define ROOT_PATH al inicio de la macro).");

setBatchMode(true);
print("Buscando stacks *_clean_quant.tif en " + raiz + " y sus subcarpetas...");
print("Carpetas que se saltan (por nombre): " + EXCLUDE_DIRS);
resultado = buscarStacks(raiz);
lineas = split(resultado, "\n");

lista = newArray(0);
extNom = newArray(0); extCnt = newArray(0);
nCarpetas = 0; nExcl = 0; nArchivos = 0; nFuera = 0;
for (q = 0; q < lineas.length; q++) {
    campos = split(lineas[q], "\t");
    tipo = campos[0];
    if (tipo == "Q") {
        lista = Array.concat(lista, campos[1]);
    } else if (tipo == "N") {
        nFuera++;
        if (VERBOSE) print("   (ignorado: no esta en una carpeta '" + INPUT_SUBDIR + "') " + campos[1]);
    } else if (tipo == "X") {
        nExcl++;
        if (VERBOSE) print("   (excluida) " + campos[1]);
    } else if (tipo == "D") {
        nCarpetas++;
        if (VERBOSE) print("   [carpeta] " + campos[1] + "  ->  " + campos[2] + " archivos, " + campos[3] + " subcarpetas, " + campos[4] + " stacks");
    } else if (tipo == "E") {
        nArchivos++;
        hallada = false;
        for (k = 0; k < extNom.length; k++) {
            if (extNom[k] == campos[1]) { extCnt[k] = extCnt[k] + 1; hallada = true; }
        }
        if (!hallada) { extNom = Array.concat(extNom, campos[1]); extCnt = Array.concat(extCnt, 1); }
    }
}
print("Carpetas visitadas: " + nCarpetas + " | excluidas: " + nExcl + " | archivos vistos: " + nArchivos + " | stacks fuera de '" + INPUT_SUBDIR + "': " + nFuera);
if (extNom.length > 0) {
    resumenExt = "";
    for (k = 0; k < extNom.length && k < 10; k++) resumenExt = resumenExt + extNom[k] + ": " + extCnt[k] + "   ";
    print("Extensiones vistas: " + resumenExt);
}
print("Stacks *_clean_quant.tif encontrados en '" + INPUT_SUBDIR + "': " + lista.length + " (fuente de ROI: " + ROI_SOURCE + ", medir: " + DO_MEASURE + ")");
if (lista.length == 0) {
    print("\nNo se encontro ningun stack *_clean_quant.tif. Posibles causas:");
    print("  - Todavia no se ha ejecutado el paso 1 (01_limpieza_lote.ijm), que deja los stacks en '" + INPUT_SUBDIR + "'.");
    print("  - Los stacks estan en otra carpeta: cambia INPUT_SUBDIR (\"\" = cualquiera). Stacks fuera de ella: " + nFuera + ".");
    print("  - La carpeta elegida no es la que crees (mira 'Carpeta recibida' arriba).");
    print("  - Los stacks estan solo dentro de carpetas excluidas (" + EXCLUDE_DIRS + ").");
    print("  - Los nombres no terminan en _clean_quant.tif (mira 'Extensiones vistas').");
    raizLista = getFileList(raiz);
    print("Contenido directo de la carpeta (" + raizLista.length + " entradas):");
    for (k = 0; k < raizLista.length && k < 20; k++) print("    " + raizLista[k]);
    exit("No se encontro ningun stack *_clean_quant.tif en:\n" + raiz + "\n\nRevisa la ventana Log (Window > Log): ahi se explica que vio la macro.");
}

// ---- StarDist: comprobaciones previas (solo si puede hacer falta) ----
testigo = raiz + "stardist_estado.txt";
hayStarDist = false;
if (ROI_SOURCE == "auto" || ROI_SOURCE == "stardist") {
    List.setCommands;
    if (List.get("StarDist 2D") == "") {
        print("StarDist: el plugin NO esta instalado en este Fiji (no aparece el comando 'StarDist 2D').");
        print("   Los stacks que tengan <nombre>_rois_labels.tif se procesan igual; los demas se omiten.");
    } else {
        // OJO: el operador && de las macros NO hace cortocircuito (evalua siempre los dos lados), asi que
        // File.openAsString solo puede llamarse dentro de un if anidado que ya sepa que el archivo existe.
        falloPrevio = false;
        if (File.exists(testigo)) {
            if (indexOf(File.openAsString(testigo), "INICIADO") >= 0) falloPrevio = true;
        }
        if (falloPrevio) {
            print("StarDist: una ejecucion anterior se detuvo al llamar al plugin; se omite para no repetir el fallo.");
            print("   Revisa Window > Console. Cuando lo hayas resuelto, borra " + testigo);
        } else {
            hayStarDist = true;
            print("StarDist: plugin instalado (esta macro no puede comprobar que funcione hasta usarlo).");
        }
    }
}

logPath = raiz + "rois_lote_log.csv";
if (File.exists(logPath)) logTxt = File.openAsString(logPath);
else logTxt = "fecha,archivo,estado,objetos_detectados,rois,descartados_por_area,area_minima_px,frames_medidos,segundos\n";

nOk = 0; nOmit = 0;
tLote = getTime();
for (k = 0; k < lista.length; k++) {
    path = lista[k];
    rel = substring(path, lengthOf(raiz));
    print("\n[" + (k + 1) + "/" + lista.length + "] " + rel);
    if (ultimaFilaInterrumpida(logTxt, rel)) {
        nOmit++;
        print("   omitido: la ejecucion anterior se detuvo con este archivo (puede estar danado).");
        print("   Para reintentarlo, borra su fila INTERRUMPIDO de " + File.getName(logPath) + ".");
        continue;
    }
    File.saveString(logTxt + ahora() + ",\"" + rel + "\",\"INTERRUMPIDO: la macro se detuvo con este archivo\",,,,,,\n", logPath);
    estado = procesarUno(path, hayStarDist, testigo);
    run("Close All");
    roiManager("reset");
    call("java.lang.System.gc");
    fecha = ahora();
    if (startsWith(estado, "OK|")) {
        c = split(estado, "|");     // OK, objetos, rois, descartados, area_minima, frames, segundos
        nOk++;
        if (c[7] == "existente") print("   OK | " + c[2] + " ROI (RoiSet existente) | " + c[5] + " frames medidos | " + c[6] + " s");
        else print("   OK | " + c[2] + " ROI de " + c[1] + " objetos, fuente " + c[7] + " (" + c[3] + " descartados por area < " + c[4] + " px) | " + c[5] + " frames medidos | " + c[6] + " s");
        logTxt = logTxt + fecha + ",\"" + rel + "\",OK," + c[1] + "," + c[2] + "," + c[3] + "," + c[4] + "," + c[5] + "," + c[6] + "\n";
    } else {
        nOmit++;
        print("   " + estado);
        if (!startsWith(estado, "omitido: ya existe"))
            logTxt = logTxt + fecha + ",\"" + rel + "\",\"" + estado + "\",,,,,,\n";     // entre comillas: el motivo puede llevar comas
    }
    File.saveString(logTxt, logPath);          // se guarda tras cada stack: permite retomar
}
print("\nLote terminado en " + d2s((getTime() - tLote) / 1000, 0) + " s: " + nOk + " procesados, " + nOmit + " omitidos.");
print("Log: " + logPath);
setBatchMode(false);

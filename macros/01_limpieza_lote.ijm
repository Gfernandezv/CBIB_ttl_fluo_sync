// ============================================================================
// 01_limpieza_lote.ijm  -  Limpieza de TODOS los stacks TIFF de una carpeta
//                          y sus subcarpetas (ImageJ/Fiji)
//
// Aplica a cada stack el mismo procesamiento que 01_limpieza_stack.ijm
// (validado contra Python):
//   1. Picos brillantes : mediana 3x3; se reemplaza si supera la mediana en mas
//                         de HOT_K x sigma robusto, calculado en CADA frame.
//   2. Deriva           : correlacion cruzada (FD Math) por bloques de BIN_SIZE
//                         frames, desplazamientos ENTEROS, mediana de la traza.
//   3. Alinear y recortar a la region comun (sin interpolar).
//   4. Vista            : fondo plano + suavizado (solo para mirar).
//
// Por cada TIFF crea, en <carpeta del TIFF>/<OUT_SUBDIR>/ :
//   <nombre>_clean_quant.tif        stack alineado, intensidades conservadas
//   <nombre>_clean_view.tif         stack para mirar (32 bits); no medir en el
//   <nombre>_drift_shifts_ij.csv    dy, dx por frame y marca de transicion
// y en la carpeta raiz: limpieza_lote_log.csv, un DIARIO con fecha: cada ejecucion
// anade una fila por TIFF evaluado (no se anaden las filas de "ya existe").
//
// QUE SE SALTA (y por que):
//   - Archivos con "_clean_" en el nombre: son salidas de esta macro.
//   - Carpetas de EXCLUDE_DIRS (por nombre): p. ej. "ROIs", que en el proyecto
//     contiene recortes de 20x20 px que tambien son stacks, y la propia carpeta
//     de salida.
//   - Imagenes que no son stacks o son muy pequenas (MIN_FRAMES, MIN_SIZE_PX).
//   - Imagenes que no son de 16 bits (solo se valido con 16 bits).
//   - Stacks que ya tienen resultados (a menos que OVERWRITE = true): se puede
//     interrumpir y volver a ejecutar; retoma donde quedo.
//
// Uso interactivo: ejecutar y elegir la carpeta raiz. Si el selector de carpetas da problemas,
//                  escribe la ruta en ROOT_PATH (arriba) y se omite el selector.
// Diagnostico: si no encuentra TIFF, el Log (Window > Log) muestra la carpeta exacta que recibio,
//              cada carpeta visitada, las extensiones que vio y la causa mas probable.
// Uso por lotes  : argumento "raiz|modo_deriva|guardar_vista" (los dos ultimos
//                  opcionales; modo = estimar | ninguna; vista = 1 | 0).
//
// VALIDADO en ImageJ 1.54g (headless) sobre dos stacks reales y un arbol de carpetas con
// casos de omision: la salida es identica a la de 01_limpieza_stack.ijm y, frente a
// Python, picos identicos pixel a pixel y deriva con diferencias de 1 px en ~10 % de los
// frames. Memoria: un stack de 305 MB corrio con 1 GB de memoria de Java.
//
// DISENO (importante): ninguna funcion de esta macro escribe variables globales. Los resultados se
// devuelven como texto y los arreglos se pasan como parametro. Motivo: segun el camino por el que se
// ejecute la macro (se comprobo que ocurre con algunos lanzadores de Fiji), las declaraciones "var" no se
// registran como globales; entonces una asignacion dentro de una funcion, como  lista = Array.concat(...),
// crea una variable LOCAL y la global queda vacia, mientras que los incrementos (n++) si funcionan. La
// version anterior de esta macro dependia de ello y fallaba en silencio ("TIFF encontrados: 0").
//
// ROBUSTEZ: una macro no puede capturar errores. Un archivo vacio se omite con un aviso. Un archivo
// corrupto, o un error grave (p. ej. falta de memoria), detiene la macro; pero antes de procesar cada
// archivo se deja en el log una fila provisional INTERRUMPIDO. Al volver a ejecutar, ese archivo se
// salta (y se explica por que) y se continua con el resto, en vez de bloquearse siempre en el mismo.
// Para reintentarlo, borra su fila INTERRUMPIDO del log.
// ============================================================================

// ---------- Parametros ----------
var HOT_K        = 5.0;     // umbral de picos, en sigmas robustos
var N_REF        = 20;      // frames iniciales usados como referencia (sin movimiento)
var BIN_SIZE     = 2;       // frames por bloque al estimar la deriva
var MEDIAN_BINS  = 3;       // ancho del filtro de mediana sobre la traza
var MIN_SHIFT_PX = 2;       // solo se corrige si la deriva maxima llega a este valor
var SIGMA_BG     = 60;      // sigma (px) del fondo gaussiano (solo vista)
var SIGMA_XY     = 1;       // suavizado espacial (px)  (solo vista)
var SIGMA_T      = 1;       // suavizado temporal (frames) (solo vista)
var DRIFT_MODE   = "estimar";   // "estimar" | "ninguna"
var SAVE_VIEW    = true;    // la vista pesa ~2 veces el stack cuantitativo
var OUT_SUBDIR   = "limpieza_ij";
var EXCLUDE_DIRS = "ROIs,limpieza_ij,procesado,batch_analysis";   // nombres separados por coma
var MIN_FRAMES   = 30;
var MIN_SIZE_PX  = 100;
var ONLY_16BIT   = true;
var OVERWRITE    = false;
var ROOT_PATH    = "";      // ruta de la carpeta raiz; si se deja vacia se abre el selector de carpetas.
                            // Se puede pegar tal cual, con barras invertidas: C:\Users\...\Proc_data
var VERBOSE      = true;    // imprime en el Log cada carpeta visitada y por que se descarta cada archivo

// Ninguna funcion de esta macro escribe variables globales: los resultados se DEVUELVEN como texto
// y los arreglos se pasan como parametro. (Solo la prueba interna de mas abajo usa globales, a proposito.)
var _gArr, _gNum, _gTxt;            // prueba interna: ver "probarGlobales"

// ---------- Funciones auxiliares ----------
function med3(a, b, c) {
    if (a > b) { t = a; a = b; b = t; }
    if (b > c) { t = b; b = c; c = t; }
    if (a > b) { t = a; a = b; b = t; }
    return b;
}

// mediana exacta a partir de un histograma de ancho de bin 1 que empieza en lo
// (mismo criterio que numpy.median: promedio de los dos centrales si n es par)
function medianaHist(cnt, lo) {
    n = 0;
    for (i = 0; i < cnt.length; i++) n += cnt[i];
    if (n % 2 == 1) { r1 = (n - 1) / 2; r2 = r1; } else { r1 = n / 2 - 1; r2 = n / 2; }
    acum = 0; v1 = lo; v2 = lo; tengo1 = false;
    for (i = 0; i < cnt.length; i++) {
        acum += cnt[i];
        if (!tengo1 && acum > r1) { v1 = lo + i; tengo1 = true; }
        if (acum > r2) { v2 = lo + i; break; }
    }
    return (v1 + v2) / 2;
}

// sigma robusto (1,4826 x MAD) del residuo pixel - mediana 3x3 de UN frame.
// El residuo de una imagen entera es entero, asi que los histogramas son exactos.
function sigmaFrame(stackTitle, idx) {
    selectImage(stackTitle);
    setSlice(idx);
    run("Duplicate...", "title=f_res");                  // solo el frame actual
    run("32-bit");
    run("Duplicate...", "title=m_res");
    run("Median...", "radius=1");
    imageCalculator("Subtract create 32-bit", "f_res", "m_res");
    rename("r_res");
    getStatistics(area, mean, mn, mx);
    lo = floor(mn); nb = floor(mx) - lo + 1;
    getHistogram(vals, cnt, nb, lo, lo + nb);
    med = medianaHist(cnt, lo);
    run("Subtract...", "value=" + d2s(med, 4));          // |residuo - mediana|, en unidades de 1/2 tras duplicar
    run("Abs");
    run("Multiply...", "value=2");                       // asi los valores x.5 pasan a enteros
    getStatistics(area, mean, mn2, mx2);
    nb2 = floor(mx2) + 1;
    getHistogram(vals2, cnt2, nb2, 0, nb2);
    mad = medianaHist(cnt2, 0) / 2;
    close("f_res"); close("m_res"); close("r_res");
    return 1.4826 * mad;
}

function pasabanda(src, dst) {
    selectImage(src); run("Duplicate...", "title=g_bajo"); run("32-bit"); run("Gaussian Blur...", "sigma=2");
    selectImage(src); run("Duplicate...", "title=g_alto"); run("32-bit"); run("Gaussian Blur...", "sigma=20");
    imageCalculator("Subtract create 32-bit", "g_bajo", "g_alto");
    rename(dst);
    close("g_bajo"); close("g_alto");
    selectImage(dst);
    run("Canvas Size...", "width=1024 height=1024 position=Top-Left zero");
}

function estimarDeriva(stackTitle, n, ody, odx) {
    selectImage(stackTitle);
    run("Z Project...", "start=1 stop=" + N_REF + " projection=[Sum Slices]");   // 32 bits (Average daria 16 bits redondeados)
    rename("ref_sum");
    pasabanda("ref_sum", "ref_bp");
    close("ref_sum");

    nb = floor((n + BIN_SIZE - 1) / BIN_SIZE);
    tdy = newArray(nb); tdx = newArray(nb);
    for (b = 0; b < nb; b++) {
        s = b * BIN_SIZE + 1;
        e = minOf(s + BIN_SIZE - 1, n);
        selectImage(stackTitle);
        run("Z Project...", "start=" + s + " stop=" + e + " projection=[Sum Slices]");
        rename("mv_sum");
        pasabanda("mv_sum", "mv_bp");
        run("FD Math...", "image1=ref_bp operation=Correlate image2=mv_bp result=cc do");
        getStatistics(area, mean, min, max);
        setThreshold(max, max);
        run("Create Selection");
        getSelectionBounds(sx, sy, sw, sh);
        tdy[b] = sy - 512;
        tdx[b] = sx - 512;
        close("cc"); close("mv_bp"); close("mv_sum");
    }
    close("ref_bp");

    fdy = newArray(nb); fdx = newArray(nb);
    for (b = 0; b < nb; b++) {
        if (MEDIAN_BINS > 1) {
            lo = maxOf(b - 1, 0); hi = minOf(b + 1, nb - 1);
            fdy[b] = med3(tdy[lo], tdy[b], tdy[hi]);
            fdx[b] = med3(tdx[lo], tdx[b], tdx[hi]);
        } else { fdy[b] = tdy[b]; fdx[b] = tdx[b]; }
    }
    for (f = 0; f < n; f++) {
        b = minOf(floor(f / BIN_SIZE), nb - 1);
        ody[f] = fdy[b]; odx[f] = fdx[b];
    }
}

function marcarTransiciones(dy, dx, flag, n) {
    for (j = 0; j < n - 1; j++) {
        salto = maxOf(abs(dy[j + 1] - dy[j]), abs(dx[j + 1] - dx[j]));
        if (salto > 3)
            for (q = maxOf(0, j - 1); q <= minOf(n - 1, j + 2); q++) flag[q] = 1;
    }
}

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

// ---------- Prueba interna ----------
// Informa si este Fiji actualiza variables globales (numero, texto, arreglo) desde una funcion.
// Es solo diagnostico: la macro no depende de ello (0 aqui no es un error).
function llenarGlobales() {
    _gNum = 5;
    _gTxt = "hola";
    _gArr = newArray(3);
    _gArr[1] = 7;
}

function probarGlobales() {
    _gNum = 0; _gTxt = ""; _gArr = newArray(0);
    llenarGlobales();
    return "numero: " + (_gNum == 5) + " | texto: " + (_gTxt == "hola") + " | arreglo: " + (_gArr.length == 3) + "   (1 = funciona, 0 = no)";
}

// ---------- Busqueda recursiva ----------
function carpetaExcluida(nombre) {
    // nombre llega con "/" final; se compara sin el y sin distinguir mayusculas (Windows no las distingue)
    if (endsWith(nombre, "/") || endsWith(nombre, "\\")) nombre = substring(nombre, 0, lengthOf(nombre) - 1);
    ex = split(EXCLUDE_DIRS, ",");
    for (i = 0; i < ex.length; i++) if (toLowerCase(String.trim(ex[i])) == toLowerCase(nombre)) return true;
    return false;
}

function extensionDe(nm) {
    pt = lastIndexOf(nm, ".");
    if (pt < 0) return "(sin extension)";
    return toLowerCase(substring(nm, pt));
}

// Devuelve TODO lo que vio, como texto (una linea por hecho, campos separados por tabulador):
//   T <ruta>                     TIFF a procesar
//   C <ruta>                     TIFF ignorado porque contiene _clean_ (es una salida)
//   X <ruta>                     carpeta excluida
//   E <extension>                un archivo (para contar extensiones)
//   D <ruta> <arch> <sub> <tif>  carpeta visitada y sus conteos
// Solo usa variables locales y valores devueltos: no depende de variables globales.
function buscarTiffs(dir) {
    nombres = getFileList(dir);
    acum = "";
    nArch = 0; nSub = 0; nTif = 0;
    for (i = 0; i < nombres.length; i++) {
        nm = nombres[i];
        if (endsWith(nm, "/") || endsWith(nm, "\\")) {
            if (carpetaExcluida(nm)) {
                acum = acum + "X\t" + dir + nm + "\n";
            } else {
                nSub++;
                acum = acum + buscarTiffs(dir + nm);
            }
        } else {
            nArch++;
            acum = acum + "E\t" + extensionDe(nm) + "\n";
            bajo = toLowerCase(nm);
            if (endsWith(bajo, ".tif") || endsWith(bajo, ".tiff")) {
                if (indexOf(nm, "_clean_") >= 0) {
                    acum = acum + "C\t" + dir + nm + "\n";
                } else {
                    acum = acum + "T\t" + dir + nm + "\n";
                    nTif++;
                }
            }
        }
    }
    acum = acum + "D\t" + dir + "\t" + nArch + "\t" + nSub + "\t" + nTif + "\n";
    return acum;
}

// ---------- Procesar un stack ----------
// Devuelve "OK|frames|ancho|alto|ancho_sal|alto_sal|deriva_max|frames_transicion|segundos" o el motivo
// por el que se omitio. No usa exit(): un stack problematico no detiene el lote.
function limpiarUno(path) {
    t0 = getTime();
    dir = File.getDirectory(path);
    nombre = File.getName(path);
    stem = substring(nombre, 0, lastIndexOf(nombre, "."));
    outDir = dir + OUT_SUBDIR + "/";
    salidaQuant = outDir + stem + "_clean_quant.tif";
    if (!OVERWRITE && File.exists(salidaQuant)) return "omitido: ya existe " + OUT_SUBDIR + "/" + stem + "_clean_quant.tif";
    if (File.length(path) == 0) return "omitido: el archivo esta vacio (0 bytes)";

    open(path);
    rename("raw");
    getDimensions(w, h, ch, sl, fr);
    if (ch >= MIN_FRAMES && sl == 1 && fr == 1) {
        // serie temporal guardada con los frames como "canales" (p. ej. por tifffile): se trata como frames
        print("   AVISO: " + ch + " canales y ningun otro eje; se trata como una serie de " + ch + " frames");
        Stack.setDimensions(1, ch, 1);
        getDimensions(w, h, ch, sl, fr);
    }
    nFrames = nSlices;
    if (ch > 1) return "omitido: tiene " + ch + " canales (imagen multicanal)";
    if (nFrames < MIN_FRAMES) return "omitido: " + nFrames + " frame(s), no es un stack (minimo " + MIN_FRAMES + ")";
    if (w < MIN_SIZE_PX || h < MIN_SIZE_PX) return "omitido: imagen muy pequena (" + w + "x" + h + ", minimo " + MIN_SIZE_PX + ")";
    if (ONLY_16BIT && bitDepth() != 16) return "omitido: " + bitDepth() + " bits (solo se valido con 16 bits)";
    if (sl == 1 && fr > 1) Stack.setDimensions(1, fr, 1);     // tratar los frames como slices
    n = nSlices;
    File.makeDirectory(outDir);

    // 1. Picos brillantes: sigma por frame (puede variar mucho entre frames); cada frame usa su propio umbral
    for (i = 1; i <= n; i++) {
        sg = sigmaFrame("raw", i);
        selectImage("raw");
        setSlice(i);
        run("Remove Outliers...", "radius=1 threshold=" + d2s(HOT_K * sg, 4) + " which=Bright");
    }

    // 2. Deriva (los arreglos se crean aqui y las funciones los rellenan)
    dyF = newArray(n); dxF = newArray(n); flagF = newArray(n);
    aplicar = false;
    dyMin = 0; dyMax = 0; dxMin = 0; dxMax = 0; maxAbs = 0; nTr = 0;
    if (DRIFT_MODE == "estimar") {
        estimarDeriva("raw", n, dyF, dxF);
        marcarTransiciones(dyF, dxF, flagF, n);
        Array.getStatistics(dyF, dyMin, dyMax);
        Array.getStatistics(dxF, dxMin, dxMax);
        maxAbs = maxOf(maxOf(abs(dyMin), abs(dyMax)), maxOf(abs(dxMin), abs(dxMax)));
        if (maxAbs >= MIN_SHIFT_PX) aplicar = true;
        else for (i = 0; i < n; i++) flagF[i] = 0;      // sin correccion no hay transiciones que marcar
        for (i = 0; i < n; i++) nTr += flagF[i];
    }
    if (aplicar) {
        selectImage("raw");
        for (i = 1; i <= n; i++) {
            if (dxF[i - 1] != 0 || dyF[i - 1] != 0) {
                setSlice(i);
                run("Translate...", "x=" + dxF[i - 1] + " y=" + dyF[i - 1] + " interpolation=None slice");
            }
        }
        yLo = maxOf(0, dyMax); yHi = h + minOf(0, dyMin);
        xLo = maxOf(0, dxMax); xHi = w + minOf(0, dxMin);
        if (yHi - yLo < 1 || xHi - xLo < 1) return "omitido: la deriva (" + maxAbs + " px) es mayor que el campo";
        makeRectangle(xLo, yLo, xHi - xLo, yHi - yLo);
        run("Crop");
    }
    if (DRIFT_MODE == "estimar") {
        tabla = "frame,dy_px,dx_px,transition\n";
        for (i = 0; i < n; i++) tabla = tabla + i + "," + dyF[i] + "," + dxF[i] + "," + flagF[i] + "\n";
        File.saveString(tabla, outDir + stem + "_drift_shifts_ij.csv");
    }

    // 3. Stack cuantitativo
    selectImage("raw");
    setSlice(1);
    w2 = getWidth(); h2 = getHeight();
    saveAs("Tiff", salidaQuant);
    rename("quant");

    // 4. Vista (solo para mirar)
    if (SAVE_VIEW) {
        selectImage("quant");
        run("Z Project...", "projection=[Sum Slices]");
        rename("media");
        run("Divide...", "value=" + n);
        run("Gaussian Blur...", "sigma=" + SIGMA_BG);
        getStatistics(area, bgMedia);
        selectImage("quant");
        run("Duplicate...", "title=view duplicate");
        run("32-bit");
        imageCalculator("Subtract stack", "view", "media");
        run("Add...", "value=" + bgMedia + " stack");
        run("Gaussian Blur 3D...", "x=" + SIGMA_XY + " y=" + SIGMA_XY + " z=" + SIGMA_T);
        resetMinAndMax();
        saveAs("Tiff", outDir + stem + "_clean_view.tif");
    }
    seg = round((getTime() - t0) / 1000);
    return "OK|" + nFrames + "|" + w + "|" + h + "|" + w2 + "|" + h2 + "|" + maxAbs + "|" + nTr + "|" + seg;
}

// ---------- Programa principal ----------
arg = getArgument();
if (arg != "") {
    p = split(arg, "|");
    raiz = p[0];
    if (p.length > 1) DRIFT_MODE = p[1];
    if (p.length > 2) SAVE_VIEW = (p[2] == "1");
} else if (ROOT_PATH != "") {
    raiz = ROOT_PATH;
} else {
    raiz = getDirectory("Elige la carpeta raiz (se buscan TIFF en ella y en sus subcarpetas)");
}
print("Carpeta recibida: [" + raiz + "]");
if (raiz == "" || raiz == "null") exit("No se eligio ninguna carpeta (el selector devolvio un valor vacio).\nDefine ROOT_PATH al inicio de la macro con la ruta de la carpeta.");
raiz = replace(raiz, "\\", "/");                 // ImageJ acepta "/" tambien en Windows; asi las rutas son uniformes
if (!endsWith(raiz, "/")) raiz = raiz + "/";
if (!File.isDirectory(raiz)) exit("No existe la carpeta: " + raiz + "\nRevisa la ruta (o define ROOT_PATH al inicio de la macro).");
print("Prueba interna de variables globales -> " + probarGlobales());

setBatchMode(true);
print("Buscando TIFF en " + raiz + " y sus subcarpetas...");
print("Carpetas que se saltan (por nombre): " + EXCLUDE_DIRS);
resultado = buscarTiffs(raiz);
lineas = split(resultado, "\n");

lista = newArray(0);
extNom = newArray(0); extCnt = newArray(0);
nCarpetas = 0; nExcl = 0; nClean = 0; nArchivos = 0;
for (q = 0; q < lineas.length; q++) {
    campos = split(lineas[q], "\t");
    tipo = campos[0];
    if (tipo == "T") {
        lista = Array.concat(lista, campos[1]);
    } else if (tipo == "C") {
        nClean++;
        if (VERBOSE) print("   (ignorado: es una salida, contiene _clean_) " + campos[1]);
    } else if (tipo == "X") {
        nExcl++;
        if (VERBOSE) print("   (excluida) " + campos[1]);
    } else if (tipo == "D") {
        nCarpetas++;
        if (VERBOSE) print("   [carpeta] " + campos[1] + "  ->  " + campos[2] + " archivos, " + campos[3] + " subcarpetas, " + campos[4] + " TIFF");
    } else if (tipo == "E") {
        nArchivos++;
        hallada = false;
        for (k = 0; k < extNom.length; k++) {
            if (extNom[k] == campos[1]) { extCnt[k] = extCnt[k] + 1; hallada = true; }
        }
        if (!hallada) { extNom = Array.concat(extNom, campos[1]); extCnt = Array.concat(extCnt, 1); }
    }
}
print("Carpetas visitadas: " + nCarpetas + " | excluidas: " + nExcl + " | archivos vistos: " + nArchivos + " | salidas ignoradas (_clean_): " + nClean);
if (extNom.length > 0) {
    resumenExt = "";
    for (k = 0; k < extNom.length && k < 10; k++) resumenExt = resumenExt + extNom[k] + ": " + extCnt[k] + "   ";
    print("Extensiones vistas: " + resumenExt);
}
print("TIFF encontrados: " + lista.length + " (modo de deriva: " + DRIFT_MODE + ", guardar vista: " + SAVE_VIEW + ")");
if (lista.length == 0) {
    print("\nNo se encontro ningun TIFF. Posibles causas:");
    print("  - La carpeta elegida no es la que crees (mira 'Carpeta recibida' arriba).");
    print("  - Los TIFF estan solo dentro de carpetas excluidas (" + EXCLUDE_DIRS + ").");
    print("  - Los archivos no terminan en .tif o .tiff (mira 'Extensiones vistas').");
    print("  - Los nombres contienen '_clean_' (se tratan como salidas).");
    print("  - La ruta no es accesible desde Fiji (p. ej. unidad de red o de la nube sin conectar).");
    raizLista = getFileList(raiz);
    print("Contenido directo de la carpeta (" + raizLista.length + " entradas):");
    for (k = 0; k < raizLista.length && k < 20; k++) print("    " + raizLista[k]);
    exit("No se encontro ningun TIFF en:\n" + raiz + "\n\nRevisa la ventana Log (Window > Log): ahi se explica que vio la macro.");
}

logPath = raiz + "limpieza_lote_log.csv";
if (File.exists(logPath)) logTxt = File.openAsString(logPath);
else logTxt = "fecha,archivo,estado,frames,entrada,salida,deriva_max_px,frames_transicion,segundos\n";

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
    estado = limpiarUno(path);
    run("Close All");
    call("java.lang.System.gc");
    fecha = ahora();
    if (startsWith(estado, "OK|")) {
        c = split(estado, "|");     // OK, frames, ancho, alto, ancho_sal, alto_sal, deriva, transiciones, segundos
        nOk++;
        print("   OK | " + c[1] + " frames | " + c[2] + "x" + c[3] + " -> " + c[4] + "x" + c[5] + " | deriva max " + c[6] + " px | " + c[7] + " frames en transicion | " + c[8] + " s");
        logTxt = logTxt + fecha + ",\"" + rel + "\",OK," + c[1] + "," + c[2] + "x" + c[3] + "," + c[4] + "x" + c[5] + "," + c[6] + "," + c[7] + "," + c[8] + "\n";
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

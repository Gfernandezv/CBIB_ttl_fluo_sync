// ============================================================================
// 01_limpieza_stack.ijm  -  Limpieza de un stack de fluorescencia (ImageJ/Fiji)
//
// Replica el flujo de Python (stack_cleanup.py), etapa por etapa:
//   1. Picos brillantes : mediana 3x3 por pixel; se reemplaza si supera la
//                         mediana en mas de HOT_K x sigma robusto, calculado
//                         en cada frame.
//   2. Deriva           : correlacion cruzada (FD Math) entre cada bloque de
//                         BIN_SIZE frames y una referencia; desplazamientos
//                         ENTEROS; mediana de la traza.
//   3. Alinear y recortar a la region comun (sin interpolar).
//   4. Vista            : fondo plano + suavizado (solo para mirar).
//
// Salidas (en la carpeta elegida):
//   <nombre>_clean_quant.tif        stack alineado, intensidades conservadas
//   <nombre>_clean_view.tif         stack para mirar (32 bits); no medir en el
//   <nombre>_drift_shifts_ij.csv    dy, dx por frame y marca de transicion
//
// DISENO (importante): ninguna funcion de esta macro escribe variables globales. Los resultados se
// devuelven como texto y los arreglos se pasan como parametro. Motivo: segun el camino por el que se
// ejecute la macro (se comprobo que ocurre con algunos lanzadores de Fiji), las declaraciones "var" no se
// registran como globales; entonces una asignacion dentro de una funcion, como  lista = Array.concat(...),
// crea una variable LOCAL y la global queda vacia, mientras que los incrementos (n++) si funcionan. La
// version anterior de esta macro dependia de ello y fallaba en silencio (se detenia tras "2/4 Deriva").
//
// Uso interactivo: ejecutar y elegir el TIFF crudo y la carpeta de salida.
// Uso por lotes  : argumento "entrada|salida|modo|csv" (modo y csv opcionales).
//
// VALIDADO en ImageJ 1.54g (headless) contra los resultados de Python:
//   - Picos: identico pixel a pixel en dos stacks reales: uno con sigma de ruido
//     constante y otro donde varia de 3 a 9 entre frames (por eso se calcula por frame).
//   - Deriva: ~90 % de los frames identicos; el resto difiere en 1 px.
//     Con DRIFT_MODE = "csv" (desplazamientos de Python) el stack cuantitativo
//     es identico pixel a pixel al de Python.
// ============================================================================

// ---------- Parametros ----------
HOT_K        = 5.0;      // umbral de picos, en sigmas robustos
N_REF        = 20;       // frames iniciales usados como referencia (sin movimiento)
BIN_SIZE     = 2;        // frames por bloque al estimar la deriva
MEDIAN_BINS  = 3;        // ancho del filtro de mediana sobre la traza
MIN_SHIFT_PX = 2;        // solo se corrige si la deriva maxima llega a este valor
SIGMA_BG     = 60;       // sigma (px) del fondo gaussiano (solo vista)
SIGMA_XY     = 1;        // suavizado espacial (px)  (solo vista)
SIGMA_T      = 1;        // suavizado temporal (frames) (solo vista)
DRIFT_MODE   = "estimar";// "estimar" | "csv" | "ninguna"
DRIFT_CSV    = "";       // CSV de Python (frame,dy_px,dx_px,transition) si modo = "csv"
SAVE_VIEW    = true;

// Nota: ninguna funcion escribe variables globales de tipo arreglo. Los arreglos se crean en el
// programa principal y se pasan como parametro (se rellenan por referencia).

// ---------- Funciones ----------
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

// pasa-banda (Gauss sigma 2 - Gauss sigma 20), 32 bits, relleno con ceros a 1024x1024
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
        // El maximo del mapa de FD Math (origen en el centro, 512,512) da el desplazamiento
        // que hay que aplicar al bloque para alinearlo con la referencia.
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

    // mediana de la traza (ventana de 3, bordes replicados)
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

function leerDerivaCSV(ruta, n, ody, odx) {
    lineas = split(File.openAsString(ruta), "\n");
    if (lineas.length - 1 != n)
        exit("El CSV tiene " + (lineas.length - 1) + " filas y el stack " + n + " frames.");
    for (i = 0; i < n; i++) {
        c = split(lineas[i + 1], ",");      // frame, dy_px, dx_px, transition
        ody[i] = parseInt(c[1]); odx[i] = parseInt(c[2]);
    }
}

function marcarTransiciones(dy, dx, flag, n) {
    for (j = 0; j < n - 1; j++) {
        salto = maxOf(abs(dy[j + 1] - dy[j]), abs(dx[j + 1] - dx[j]));
        if (salto > 3)
            for (q = maxOf(0, j - 1); q <= minOf(n - 1, j + 2); q++) flag[q] = 1;
    }
}

// ---------- Entrada ----------
arg = getArgument();
if (arg != "") {
    p = split(arg, "|");
    inPath = p[0];
    outDir = p[1];
    if (p.length > 2) DRIFT_MODE = p[2];
    if (p.length > 3) DRIFT_CSV = p[3];
} else {
    inPath = File.openDialog("Elige el stack crudo (.tif)");
    outDir = getDirectory("Elige la carpeta de salida");
}
if (!endsWith(outDir, File.separator)) outDir = outDir + File.separator;
File.makeDirectory(outDir);
nombre = File.getName(inPath);
stem = substring(nombre, 0, lastIndexOf(nombre, "."));

setBatchMode(true);
t0 = getTime();
open(inPath);
rename("raw");
getDimensions(w, h, ch, sl, fr);
if (ch > 1) exit("Se espera un solo canal; el archivo tiene " + ch + ".");
if (bitDepth() != 16) print("AVISO: se probo con stacks de 16 bits; este es de " + bitDepth() + " bits.");
if (sl == 1 && fr > 1) Stack.setDimensions(1, fr, 1);     // tratar los frames como slices
n = nSlices;
print("Stack: " + w + " x " + h + " px, " + n + " frames, " + bitDepth() + " bits");

// ---------- 1. Picos brillantes ----------
// El sigma del ruido puede cambiar mucho de un frame a otro (en un stack real fue de 3 a 9),
// asi que se calcula en CADA frame, como en Python, y se aplica su propio umbral.
sMin = 1e9; sMax = 0;
for (i = 1; i <= n; i++) {
    sg = sigmaFrame("raw", i);
    sMin = minOf(sMin, sg); sMax = maxOf(sMax, sg);
    selectImage("raw");
    setSlice(i);
    run("Remove Outliers...", "radius=1 threshold=" + d2s(HOT_K * sg, 4) + " which=Bright");   // sin "stack": solo este frame
}
print("1/4 Picos: sigma robusto por frame entre " + d2s(sMin, 2) + " y " + d2s(sMax, 2) + " -> umbral entre " + d2s(HOT_K * sMin, 1) + " y " + d2s(HOT_K * sMax, 1) + " cuentas");

// ---------- 2. Deriva ----------
print("2/4 Deriva (modo: " + DRIFT_MODE + ")");
aplicar = false;
dyF = newArray(n); dxF = newArray(n); flagF = newArray(n);     // desplazamiento por frame y marca de transicion
if (DRIFT_MODE == "estimar") { estimarDeriva("raw", n, dyF, dxF); aplicar = true; }
if (DRIFT_MODE == "csv")     { leerDerivaCSV(DRIFT_CSV, n, dyF, dxF); aplicar = true; }
if (aplicar) {
    marcarTransiciones(dyF, dxF, flagF, n);
    Array.getStatistics(dyF, dyMin, dyMax);
    Array.getStatistics(dxF, dxMin, dxMax);
    maxAbs = maxOf(maxOf(abs(dyMin), abs(dyMax)), maxOf(abs(dxMin), abs(dxMax)));
    nTr = 0; for (i = 0; i < n; i++) nTr += flagF[i];
    print("   deriva maxima: " + maxAbs + " px | frames en transicion: " + nTr);
    if (maxAbs < MIN_SHIFT_PX) { aplicar = false; print("   menor que " + MIN_SHIFT_PX + " px: no se corrige"); }
}

if (aplicar) {
    // (+dx, +dy) mueve el contenido a la derecha/abajo, igual que scipy.ndimage.shift
    selectImage("raw");
    for (i = 1; i <= n; i++) {
        if (dxF[i - 1] != 0 || dyF[i - 1] != 0) {
            setSlice(i);
            run("Translate...", "x=" + dxF[i - 1] + " y=" + dyF[i - 1] + " interpolation=None slice");
        }
    }
    yLo = maxOf(0, dyMax); yHi = h + minOf(0, dyMin);
    xLo = maxOf(0, dxMax); xHi = w + minOf(0, dxMin);
    if (yHi - yLo < 1 || xHi - xLo < 1) exit("La deriva es mayor que el campo: no hay region comun.");
    makeRectangle(xLo, yLo, xHi - xLo, yHi - yLo);
    run("Crop");
    print("   recorte: x " + xLo + ".." + xHi + ", y " + yLo + ".." + yHi);
}

// tabla de desplazamientos (mismo formato que el CSV de Python; frame desde 0)
if (aplicar) {
    s = "frame,dy_px,dx_px,transition\n";
    for (i = 0; i < n; i++) s = s + i + "," + dyF[i] + "," + dxF[i] + "," + flagF[i] + "\n";
    File.saveString(s, outDir + stem + "_drift_shifts_ij.csv");
}

// ---------- 3. Stack cuantitativo ----------
selectImage("raw");
setSlice(1);
saveAs("Tiff", outDir + stem + "_clean_quant.tif");
rename("quant");
print("3/4 Guardado: " + stem + "_clean_quant.tif (" + getWidth() + " x " + getHeight() + " x " + nSlices + ")");

// ---------- 4. Vista (solo para mirar) ----------
if (SAVE_VIEW) {
    selectImage("quant");
    run("Z Project...", "projection=[Sum Slices]");        // 32 bits
    rename("media");
    run("Divide...", "value=" + n);
    run("Gaussian Blur...", "sigma=" + SIGMA_BG);
    getStatistics(area, bgMedia);
    selectImage("quant");
    run("Duplicate...", "title=view duplicate");
    run("32-bit");
    imageCalculator("Subtract stack", "view", "media");    // resta el fondo a cada frame
    run("Add...", "value=" + bgMedia + " stack");           // repone el nivel medio
    run("Gaussian Blur 3D...", "x=" + SIGMA_XY + " y=" + SIGMA_XY + " z=" + SIGMA_T);
    resetMinAndMax();
    saveAs("Tiff", outDir + stem + "_clean_view.tif");
    print("4/4 Guardado: " + stem + "_clean_view.tif (solo para mirar, no usar para medir)");
}

print("Listo en " + d2s((getTime() - t0) / 1000, 0) + " s. Salida: " + outDir);
setBatchMode(false);

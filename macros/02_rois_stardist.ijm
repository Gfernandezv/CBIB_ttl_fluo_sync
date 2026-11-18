// ============================================================================
// 02_rois_stardist.ijm  -  Definir ROI con StarDist 2D y guardar el RoiSet
//
// Replica stardist_rois.py:
//   1. Resumen temporal 2D del stack "para mirar" (media por defecto).
//   2. Segmentacion -> imagen de etiquetas (un numero por objeto).
//        ROI_SOURCE = "stardist"  : plugin StarDist 2D de Fiji (ver nota).
//        ROI_SOURCE = "etiquetas" : abre un TIFF de etiquetas hecho en otro
//                                   lado (p. ej. el de Python).
//   3. Etiquetas -> ROI: descarta objetos con area < MIN_AREA (automatica si es 0), deja la
//      componente mayor de cada uno, traza el contorno por los bordes de
//      pixel (el ROI reproduce la mascara exactamente) y ordena de arriba
//      hacia abajo (por el centroide; a igual y, de izquierda a derecha).
//   4. Guarda RoiSet.zip. La posicion k en el ROI Manager es ROIk.
//
// CONECTIVIDAD: para quedarse con la componente mayor de cada objeto, Analyze Particles usa conectividad
// 8 (un pixel unido solo en diagonal cuenta como parte del objeto). El modulo de Python (stardist_rois.py)
// usa conectividad 4 y lo descarta. En un stack real, 1 de 278 ROI difirio en 1 pixel por este motivo
// (el resto fue identico); si ocurre, ese ROI tiene 1 pixel mas en ImageJ y su fila de ROI_Results difiere.
//
// Los ROI estan en las coordenadas del stack ALINEADO y RECORTADO: abrelos
// sobre *_clean_quant.tif, no sobre el stack crudo.
//
// Uso interactivo: ejecutar y elegir el stack de vista (o el TIFF de etiquetas).
// Uso por lotes  : argumento "vista|salida|fuente|etiquetas" (los dos ultimos opcionales).
//
// VALIDADO en ImageJ 1.54g (headless): la conversion etiquetas -> ROI -> RoiSet
// (descarte por area, componente mayor, orden y guardado) contra Python.
// NO VALIDADO: la llamada al plugin StarDist de Fiji (ROI_SOURCE = "stardist"),
// que requiere los sitios de actualizacion CSBDeep y StarDist, y que podria dar
// una segmentacion algo distinta de la de Python.
// ============================================================================

// ---------- Parametros ----------
ROI_SOURCE  = "stardist";   // "stardist" | "etiquetas"
LABELS_PATH = "";           // TIFF de etiquetas si ROI_SOURCE = "etiquetas"
PROJECTION  = "mean";       // "mean" | "max"
MIN_AREA    = 0;            // area minima de un ROI (px). 0 = automatico: max(20, 25 % de la mediana de las areas).
                            // Un valor fijo no sirve entre stacks: 150 px dejaba 3 ROI de 278 en un stack de celulas pequenas.
// Parametros de StarDist (los del modelo 2D_versatile_fluo; coinciden con Python)
STARDIST_MODEL = "Versatile (fluorescent nuclei)";
PROB_THRESH    = 0.479071;
NMS_THRESH     = 0.3;
EXCLUDE_BOUNDARY = 0;       // 0 = no excluir objetos que tocan el borde (como en Python)

// ---------- Entrada ----------
arg = getArgument();
if (arg != "") {
    p = split(arg, "|");
    inPath = p[0];
    outDir = p[1];
    if (p.length > 2) ROI_SOURCE = p[2];
    if (p.length > 3) LABELS_PATH = p[3];
} else {
    if (ROI_SOURCE == "etiquetas") inPath = File.openDialog("Elige el TIFF de etiquetas");
    else inPath = File.openDialog("Elige el stack de vista (*_clean_view.tif)");
    outDir = getDirectory("Elige la carpeta de salida");
}
if (!endsWith(outDir, File.separator)) outDir = outDir + File.separator;
File.makeDirectory(outDir);
if (ROI_SOURCE == "etiquetas" && LABELS_PATH == "") LABELS_PATH = inPath;

setBatchMode(true);
t0 = getTime();
roiManager("reset");

// ---------- 1-2. Etiquetas ----------
if (ROI_SOURCE == "etiquetas") {
    open(LABELS_PATH);
    rename("labels");
    print("Etiquetas leidas de " + File.getName(LABELS_PATH));
} else {
    open(inPath);
    rename("vista");
    getDimensions(w, h, ch, sl, fr);
    if (sl == 1 && fr > 1) Stack.setDimensions(1, fr, 1);
    n = nSlices;
    if (PROJECTION == "max") {
        run("Z Project...", "projection=[Max Intensity]");
    } else {
        run("Z Project...", "projection=[Sum Slices]");      // 32 bits
        run("Divide...", "value=" + n);
    }
    rename("media2d");
    close("vista");
    // Llamada al plugin StarDist 2D (Plugins > StarDist > StarDist 2D). Parametros
    // tal como los registra la grabadora de macros de Fiji. SIN PROBAR en esta sesion.
    run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], "
        + "args=['input':'media2d', 'modelChoice':'" + STARDIST_MODEL + "', "
        + "'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', "
        + "'probThresh':'" + PROB_THRESH + "', 'nmsThresh':'" + NMS_THRESH + "', "
        + "'outputType':'Label Image', 'nTiles':'1', 'excludeBoundary':'" + EXCLUDE_BOUNDARY + "', "
        + "'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', "
        + "'showProbAndDist':'false'], process=[false]");
    selectImage("Label Image");            // titulo con el que el plugin crea la imagen de etiquetas
    rename("labels");
}
selectImage("labels");
getDimensions(lw, lh, lc, ls, lf);
if (ls * lf > 1) exit("La imagen de etiquetas debe ser 2D.");
getStatistics(area, mean, min, nEtiquetas);
print("Objetos detectados: " + nEtiquetas);

// area de cada etiqueta (histograma exacto) y umbral de area minima
if (bitDepth() == 16) getHistogram(vals, cnt, 65536);
else if (bitDepth() == 8) getHistogram(vals, cnt, 256);
else exit("La imagen de etiquetas debe ser de 8 o 16 bits (tiene " + bitDepth() + ").");
if (MIN_AREA > 0) {
    minArea = MIN_AREA;
} else {
    pos = newArray(0);
    for (k = 1; k <= nEtiquetas; k++) if (cnt[k] > 0) pos = Array.concat(pos, cnt[k]);
    if (pos.length == 0) exit("La imagen de etiquetas no tiene objetos.");
    Array.sort(pos);
    nPos = pos.length;
    if (nPos % 2 == 1) mediana = pos[(nPos - 1) / 2]; else mediana = (pos[nPos / 2 - 1] + pos[nPos / 2]) / 2;
    minArea = maxOf(20, 0.25 * mediana);
}
if (MIN_AREA > 0) print("Area minima de un ROI: " + d2s(minArea, 1) + " px");
else print("Area minima de un ROI: " + d2s(minArea, 1) + " px (automatica)");

// ---------- 3. Etiquetas -> ROI ----------
run("Set Measurements...", "area centroid redirect=None decimal=3");
descartados = 0;
for (k = 1; k <= nEtiquetas; k++) {
    selectImage("labels");
    run("Select None");                            // sin esto, medir/analizar se restringe al ROI anterior
    setThreshold(k, k);
    List.setMeasurements("limit");                 // medir solo los pixeles dentro del umbral
    areaK = List.getValue("Area");
    if (areaK == 0) { resetThreshold(); continue; }                  // etiqueta vacia
    if (areaK < minArea) { descartados++; resetThreshold(); continue; }

    n0 = roiManager("count");
    run("Analyze Particles...", "size=1-Infinity pixel include add");   // un ROI por componente conexa
    resetThreshold();
    n1 = roiManager("count");
    if (n1 == n0) continue;

    // quedarse con la componente de mayor area
    mejor = n0; mejorA = -1;
    for (i = n0; i < n1; i++) {
        roiManager("select", i);
        getStatistics(a);
        if (a > mejorA) { mejorA = a; mejor = i; }
    }
    for (i = n1 - 1; i > mejor; i--) { roiManager("select", i); roiManager("delete"); }
    for (i = mejor - 1; i >= n0; i--) { roiManager("select", i); roiManager("delete"); }
    idx = roiManager("count") - 1;

    // clave de orden: centroide y, luego x (rellenada con ceros para ordenar como texto)
    roiManager("select", idx);
    List.setMeasurements;
    cx = List.getValue("X"); cy = List.getValue("Y");
    roiManager("rename", IJ.pad(round(cy * 1000), 9) + "_" + IJ.pad(round(cx * 1000), 9));
}
nRoi = roiManager("count");
print("ROI: " + nRoi + " (descartados por area < " + d2s(minArea, 1) + " px: " + descartados + ")");
if (nRoi == 0) exit("No quedo ningun ROI.");

// ---------- 4. Orden y nombres ----------
roiManager("sort");
for (i = 0; i < nRoi; i++) {                       // nombre estilo ImageJ "yyyy-xxxx" (centro del cuadro)
    roiManager("select", i);
    getSelectionBounds(bx, by, bw, bh);
    roiManager("rename", IJ.pad(floor(by + bh / 2), 4) + "-" + IJ.pad(floor(bx + bw / 2), 4));
}
roiManager("deselect");
roiManager("save", outDir + "RoiSet.zip");

// tabla de ROI (centroide y area), por si se quiere cruzar con Python
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
saveAs("Results", outDir + "rois_table_ij.csv");
print("Guardado: RoiSet.zip y rois_table_ij.csv en " + outDir);
print("Listo en " + d2s((getTime() - t0) / 1000, 1) + " s");
setBatchMode(false);

// ============================================================================
// 03_medicion.ijm  -  Medir los ROI y guardar ROI_Results.csv
//
// Equivalente a stardist_rois.measure_rois / write_roi_results. Por cada ROI y
// cada frame mide Area (pixeles), IntDen y RawIntDen, y escribe la tabla con el
// mismo formato que el Multi Measure de ImageJ: una primera columna sin nombre
// con el indice de frame (desde 1) y, por cada ROI k, las columnas Area{k},
// IntDen{k} y RawIntDen{k}, en una fila por frame.
//
// NO usa el comando Multi Measure del ROI Manager: desde una macro no permite
// elegir de forma fiable el diseno "una fila por frame" y su ejecucion es
// asincrona. Aqui la suma se calcula directamente por ROI y por frame.
//
// Se mide sobre el stack CUANTITATIVO (*_clean_quant.tif), nunca sobre la vista
// suavizada. La imagen no debe tener calibracion (1 px = 1 unidad): asi
// IntDen = RawIntDen.
//
// Uso interactivo: ejecutar y elegir el stack cuantitativo, el RoiSet y la
// carpeta de salida.
// Uso por lotes  : argumento "quant|roiset|salida".
//
// VALIDADO en ImageJ 1.54g (headless) contra la medicion de Python.
// ============================================================================

// ---------- Entrada ----------
arg = getArgument();
if (arg != "") {
    p = split(arg, "|");
    quantPath = p[0];
    roiPath = p[1];
    outDir = p[2];
} else {
    quantPath = File.openDialog("Elige el stack cuantitativo (*_clean_quant.tif)");
    roiPath = File.openDialog("Elige el RoiSet (.zip)");
    outDir = getDirectory("Elige la carpeta de salida");
}
if (!endsWith(outDir, File.separator)) outDir = outDir + File.separator;
File.makeDirectory(outDir);

setBatchMode(true);
t0 = getTime();
open(quantPath);
rename("quant");
getDimensions(w, h, ch, sl, fr);
if (sl == 1 && fr > 1) Stack.setDimensions(1, fr, 1);     // tratar los frames como slices
n = nSlices;
if (bitDepth() == 32) exit("Se espera el stack cuantitativo (enteros), no la vista de 32 bits.");
getPixelSize(unidad, pw, ph);
if (pw != 1 || ph != 1)
    print("AVISO: la imagen tiene calibracion (" + pw + " x " + ph + " " + unidad + "). Se mide en pixeles igualmente.");
print("Stack: " + w + " x " + h + " px, " + n + " frames");

roiManager("reset");
roiManager("open", roiPath);
nRoi = roiManager("count");
print("ROI cargados: " + nRoi);
if (nRoi == 0) exit("El RoiSet no tiene ROI.");

// ---------- Medicion ----------
// suma[k * n + (frame - 1)] = suma de los valores del ROI k (0-based) en ese frame
suma = newArray(nRoi * n);
area = newArray(nRoi);
selectImage("quant");
for (k = 0; k < nRoi; k++) {
    roiManager("select", k);                       // el ROI se mantiene al cambiar de slice
    for (i = 1; i <= n; i++) {
        setSlice(i);
        getRawStatistics(nPix, media);
        if (i == 1) area[k] = nPix;
        suma[k * n + (i - 1)] = round(nPix * media);   // suma exacta de enteros
    }
}
run("Select None");

// ---------- Tabla ----------
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
print("Tabla: " + nResults + " filas (frames) x " + (nRoi * 3) + " columnas de medidas");
saveAs("Results", outDir + "ROI_Results.csv");
print("Guardado: ROI_Results.csv en " + outDir);
print("Listo en " + d2s((getTime() - t0) / 1000, 0) + " s");
setBatchMode(false);

<#
actualizar_rawdata.ps1
Actualiza las columnas limpieza / Stardist / ROI analysis de la hoja "Raw Data"
mirando que archivos existen en el disco. No requiere Python (PowerShell 5.1+).

Que marca como "ok":
  limpieza      -> existe  <stem>_clean_quant.tif   (carpeta limpieza_ij)
  Stardist      -> existe  <stem>_RoiSet.zip        (carpeta rois_ij)
  ROI analysis  -> existe  <stem>_ROI_Results.csv   (carpeta rois_ij)

USO (en PowerShell):
  # 1) Solo ver que hay en el disco (no necesita tabla):
  .\actualizar_rawdata.ps1 -Raiz "D:\German Fernandez\TRPM3\Fluo4"

  # 2) Actualizar la tabla: en Google Sheets, pestana "Raw Data" ->
  #    Archivo > Descargar > Valores separados por comas (.csv) y guardalo como RawData.csv
  .\actualizar_rawdata.ps1 -Raiz "D:\German Fernandez\TRPM3\Fluo4" -Tabla "C:\ruta\RawData.csv"

Resultado: RawData_actualizado.csv (misma tabla, mismo orden de filas) y
estado_disco.csv (todo lo encontrado en el disco). El original NO se modifica.
Luego copia las columnas limpieza / Stardist / ROI analysis y pegalas en la hoja.

Si el archivo no corre por la politica de ejecucion:
  powershell -ExecutionPolicy Bypass -File .\actualizar_rawdata.ps1 -Raiz "D:\..."
#>
param(
    [Parameter(Mandatory = $true)][string]$Raiz,
    [string]$Tabla = "",
    [string]$Salida = "RawData_actualizado.csv",
    [string]$Marca = "ok",
    [string]$ColImagen = "image File",
    [string]$ColLimpieza = "limpieza",
    [string]$ColStardist = "Stardist",
    [string]$ColRoi = "ROI analysis"
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path -LiteralPath $Raiz)) {
    Write-Host "No existe la carpeta raiz: $Raiz" -ForegroundColor Red
    exit 1
}

# Nombre -> clave comparable: minusculas, sin extension ni sufijos de pipeline
function Get-Clave([string]$nombre) {
    $s = $nombre.Trim().ToLower()
    $s = $s -replace '\.(tif|tiff)$', ''
    $s = $s -replace '_(clean_quant|clean_view|roiset|roi_results|rois_table_ij|drift_shifts_ij)$', ''
    $s = $s -replace '_(stack|roi)$', ''
    return $s
}

# ---- 1. Escanear disco -------------------------------------------------
Write-Host "Escaneando $Raiz ..."
$estado = @{}   # clave -> objeto

function Get-Entrada([string]$clave, [string]$stem) {
    if (-not $estado.ContainsKey($clave)) {
        $estado[$clave] = [PSCustomObject]@{
            Clave = $clave; Stem = $stem
            Limpieza = ""; Stardist = ""; RoiAnalysis = ""; Carpeta = ""
        }
    }
    return $estado[$clave]
}

$archivos = Get-ChildItem -LiteralPath $Raiz -Recurse -File -ErrorAction SilentlyContinue |
    Where-Object { $_.Name -match '_clean_quant\.tif+$|_RoiSet\.zip$|_ROI_Results\.csv$' }

foreach ($f in $archivos) {
    $n = $f.Name
    if ($f.Length -eq 0) { continue }                      # archivo vacio = no cuenta
    $clave = Get-Clave $n
    $stem = $n -replace '_(clean_quant\.tif+|RoiSet\.zip|ROI_Results\.csv)$', ''
    $e = Get-Entrada $clave $stem
    if ($n -match '_clean_quant\.tif+$') { $e.Limpieza = $Marca }
    elseif ($n -match '_RoiSet\.zip$') { $e.Stardist = $Marca }
    elseif ($n -match '_ROI_Results\.csv$') { $e.RoiAnalysis = $Marca }
    if ($e.Carpeta -eq "") { $e.Carpeta = $f.DirectoryName }
}

$lista = $estado.Values | Sort-Object Stem
$lista | Export-Csv -Path "estado_disco.csv" -NoTypeInformation -Encoding UTF8 -Delimiter ","
Write-Host ("Stacks con algun resultado en disco: {0}  (detalle en estado_disco.csv)" -f @($lista).Count)

if ($Tabla -eq "") {
    $lista | Format-Table Stem, Limpieza, Stardist, RoiAnalysis -AutoSize
    Write-Host "Sin -Tabla: no se actualizo ninguna tabla." -ForegroundColor Yellow
    exit 0
}

# ---- 2. Actualizar la tabla -------------------------------------------
if (-not (Test-Path -LiteralPath $Tabla)) {
    Write-Host "No existe la tabla: $Tabla" -ForegroundColor Red
    exit 1
}

# Delimiter explicito: en Windows en espanol el predeterminado seria ";"
$filas = Import-Csv -LiteralPath $Tabla -Delimiter ","
if (@($filas).Count -eq 0) { Write-Host "La tabla esta vacia." -ForegroundColor Red; exit 1 }

$cols = $filas[0].PSObject.Properties.Name
foreach ($c in @($ColImagen, $ColLimpieza, $ColStardist, $ColRoi)) {
    if ($cols -notcontains $c) {
        Write-Host "Falta la columna '$c'. Columnas encontradas: $($cols -join ' | ')" -ForegroundColor Red
        exit 1
    }
}

$vistos = @{}
$cambios = 0
$sinDisco = @()
foreach ($r in $filas) {
    $nombre = [string]$r.$ColImagen
    if ($nombre.Trim() -eq "") { continue }
    $clave = Get-Clave $nombre
    $vistos[$clave] = $true
    if ($estado.ContainsKey($clave)) {
        $e = $estado[$clave]
        if ($e.Limpieza -ne "" -and $r.$ColLimpieza -ne $Marca) { $r.$ColLimpieza = $Marca; $cambios++ }
        if ($e.Stardist -ne "" -and $r.$ColStardist -ne $Marca) { $r.$ColStardist = $Marca; $cambios++ }
        if ($e.RoiAnalysis -ne "" -and $r.$ColRoi -ne $Marca) { $r.$ColRoi = $Marca; $cambios++ }
    } else {
        $sinDisco += $nombre
    }
}

$filas | Export-Csv -Path $Salida -NoTypeInformation -Encoding UTF8 -Delimiter ","

# ---- 3. Resumen --------------------------------------------------------
Write-Host ""
Write-Host ("Celdas actualizadas: {0}  -> {1}" -f $cambios, $Salida) -ForegroundColor Green

if ($sinDisco.Count -gt 0) {
    Write-Host ""
    Write-Host "Filas de la tabla SIN ningun resultado en disco:" -ForegroundColor Yellow
    $sinDisco | Sort-Object -Unique | ForEach-Object { Write-Host "  $_" }
}
$huerfanos = $lista | Where-Object { -not $vistos.ContainsKey($_.Clave) }
if (@($huerfanos).Count -gt 0) {
    Write-Host ""
    Write-Host "Resultados en disco que NO estan en la tabla:" -ForegroundColor Yellow
    $huerfanos | ForEach-Object { Write-Host "  $($_.Stem)   ($($_.Carpeta))" }
}
Write-Host ""
Write-Host "Aviso: si un mismo nombre aparece en varias filas (ej. Nov052026_14_Stack), todas reciben el mismo estado."

<#
actualizar_checklist.ps1
Un solo comando para actualizar el checklist (hoja "Raw Data") despues de correr el pipeline.
Usa verificar_pipeline.ps1 (debe estar en la MISMA carpeta que este archivo).

Pasos que hace:
  1. (opcional) ordena archivos fuera de lugar     -> -Normalizar [-Aplicar]
  2. verifica cada .tif original: limpieza / Stardist / ROI analysis
  3. actualiza la tabla Raw Data y agrega los stacks nuevos al final
  4. guarda RawData_actualizado_AAAAMMDD_HHMM.csv (el original NO se toca)
  5. copia la tabla completa al portapapeles, lista para pegar en la hoja (A1 de "Raw Data")

USO (desde la carpeta Fluo4):
  1) En Google Sheets, pestana "Raw Data": Archivo > Descargar > CSV  -> guardalo como RawData.csv en Fluo4
  2) .\actualizar_checklist.ps1
     (busca RawData*.csv en la carpeta actual; o indica -Tabla "ruta\archivo.csv")
  3) En la hoja: clic en A1 de "Raw Data" > Ctrl+V. Si algo sale mal: Ctrl+Z.

Opciones:
  -Raiz "."            carpeta Fluo4 (la que contiene 1227, 1265, Wt...)
  -Tabla "ruta.csv"    tabla exportada
  -Normalizar          muestra el plan para ordenar carpetas (no mueve nada)
  -Normalizar -Aplicar mueve los archivos
  -SinPortapapeles     no copia al portapapeles
  -NoAgregar           no agrega filas de stacks que faltan en la tabla
Si Windows bloquea el script:
  powershell -ExecutionPolicy Bypass -File .\actualizar_checklist.ps1
#>
param(
    [string]$Raiz = ".",
    [string]$Tabla = "",
    [switch]$Normalizar,
    [switch]$Aplicar,
    [switch]$SinPortapapeles,
    [switch]$NoAgregar,
    [switch]$SinROIs
)

$ErrorActionPreference = "Stop"

# ---- localizar verificar_pipeline.ps1 (misma carpeta que este script) ----
$aqui = $PSScriptRoot
if ([string]::IsNullOrEmpty($aqui)) { $aqui = (Get-Location).Path }
$verif = Join-Path $aqui "verificar_pipeline.ps1"
if (-not (Test-Path -LiteralPath $verif)) {
    Write-Host "No encuentro verificar_pipeline.ps1 junto a este script ($aqui)." -ForegroundColor Red
    exit 1
}
if ((Get-Item -LiteralPath $verif).Length -eq 0) {
    Write-Host "verificar_pipeline.ps1 esta vacio (0 bytes). Vuelve a copiarlo." -ForegroundColor Red
    exit 1
}

# ---- localizar la tabla ----
if ($Tabla -eq "") {
    $cand = Get-ChildItem -LiteralPath (Get-Location).Path -File -Filter "*.csv" -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -match 'raw\s*data' -and $_.Name -notmatch 'actualizado' } |
        Sort-Object LastWriteTime -Descending
    if (@($cand).Count -eq 0) {
        Write-Host "No encontre un RawData*.csv en $((Get-Location).Path). Exporta la hoja como CSV o usa -Tabla." -ForegroundColor Red
        exit 1
    }
    $Tabla = @($cand)[0].FullName
    Write-Host "Tabla: $Tabla"
}
if (-not (Test-Path -LiteralPath $Tabla)) {
    Write-Host "No existe la tabla: $Tabla" -ForegroundColor Red
    exit 1
}

$sello = Get-Date -Format "yyyyMMdd_HHmm"
$salidaTabla = Join-Path (Get-Location).Path ("RawData_actualizado_" + $sello + ".csv")

# ---- ejecutar la verificacion ----
$p = @{ Raiz = $Raiz; Tabla = $Tabla; SalidaTabla = $salidaTabla }
if ($Normalizar) { $p["Normalizar"] = $true }
if ($Aplicar) { $p["Aplicar"] = $true }
if ($NoAgregar) { $p["NoAgregar"] = $true }
if ($SinROIs) { $p["SinROIs"] = $true }

& $verif @p
$codigo = $LASTEXITCODE
if ($codigo -eq 1) {
    Write-Host "La verificacion no pudo completarse (ver mensajes arriba)." -ForegroundColor Red
    exit 1
}

if (-not (Test-Path -LiteralPath $salidaTabla)) {
    Write-Host "No se genero $salidaTabla." -ForegroundColor Red
    exit 1
}

# ---- portapapeles: tabla completa separada por tabuladores, sin comillas ----
if (-not $SinPortapapeles) {
    $filas = @(Import-Csv -LiteralPath $salidaTabla -Delimiter ",")
    if ($filas.Count -gt 0) {
        $cols = $filas[0].PSObject.Properties.Name
        $lineas = New-Object System.Collections.ArrayList
        [void]$lineas.Add(($cols -join "`t"))
        foreach ($f in $filas) {
            $celdas = @()
            foreach ($c in $cols) {
                $v = [string]$f.$c
                $v = $v.Replace("`t", " ").Replace("`r", " ").Replace("`n", " ")
                $celdas += $v
            }
            [void]$lineas.Add(($celdas -join "`t"))
        }
        ($lineas -join "`r`n") | Set-Clipboard
        Write-Host ""
        Write-Host ("Tabla copiada al portapapeles: {0} filas + encabezado." -f $filas.Count) -ForegroundColor Green
        Write-Host "En Google Sheets: pestana 'Raw Data' > clic en A1 > Ctrl+V  (Ctrl+Z si algo sale mal)."
    }
}

Write-Host ""
Write-Host "Archivo guardado: $salidaTabla" -ForegroundColor Green
if ($Normalizar -and -not $Aplicar) {
    Write-Host "Recuerda: -Normalizar sin -Aplicar fue solo una simulacion (ver normalizacion_log.csv)." -ForegroundColor Yellow
}

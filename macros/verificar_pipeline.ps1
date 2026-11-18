<#
verificar_pipeline.ps1
Verifica la integridad del pipeline: que cada .tif original haya pasado por
  paso 1 (limpieza)  -> limpieza_ij\<stem>_clean_quant.tif (+ _clean_view.tif, _drift_shifts_ij.csv)
  paso 2 (ROI)       -> rois_ij\<stem>_RoiSet.zip
  paso 3 (medicion)  -> rois_ij\<stem>_ROI_Results.csv
y que los archivos sean coherentes entre si. Solo LEE: no modifica ni borra nada.
Requiere PowerShell 5.1+ (viene con Windows). No necesita Python.

Comprobaciones por stack:
  - existe cada salida y no esta vacia (0 bytes)
  - el RoiSet.zip se abre y contiene ROIs (.roi)
  - n.o de ROI del RoiSet = n.o de columnas Area{k} de ROI_Results.csv
  - n.o de frames: original = clean_quant = filas de ROI_Results = filas del CSV de deriva
    (se lee de la cabecera ImageJ del TIFF: "images=N"; si no esta, esa comparacion se omite)
  - orden cronologico: clean_quant mas nuevo que el original; RoiSet/ROI_Results mas nuevos que clean_quant

Estados (columna Estado):
  OK                 todo presente y coherente
  FALTA_LIMPIEZA     no hay clean_quant (el paso 1 no se ejecuto o fallo)
  FALTA_ROIS         hay limpieza pero no RoiSet
  FALTA_MEDICION     hay RoiSet pero no ROI_Results
  ERROR              archivo vacio/corrupto o incoherente (ver columna Detalle)
  Los avisos no cambian el Estado: quedan en la columna Avisos.

USO:
  .\verificar_pipeline.ps1 -Raiz "D:\German Fernandez\TRPM3\Fluo4"
  .\verificar_pipeline.ps1 -Raiz "D:\..." -SinROIs        # si solo se hizo el paso 1
  # Ademas, actualizar la hoja Raw Data (exportada como CSV desde Google Sheets):
  .\verificar_pipeline.ps1 -Raiz "D:\..." -Tabla "C:\ruta\RawData.csv"
  # Ordenar archivos que quedaron en otra carpeta (hechos con versiones distintas del pipeline):
  .\verificar_pipeline.ps1 -Raiz "D:\..." -Normalizar            # solo MUESTRA el plan, no mueve nada
  .\verificar_pipeline.ps1 -Raiz "D:\..." -Normalizar -Aplicar   # mueve los archivos
Disposicion estandar, por muestra (la carpeta del .tif original):
  <muestra>\<stem>.tif
  <muestra>\limpieza_ij\  <stem>_clean_quant.tif, _clean_view.tif, _drift_shifts_ij.csv, _rois_labels.tif
  <muestra>\rois_ij\      <stem>_RoiSet.zip, _ROI_Results.csv, _rois_table_ij.csv
-Normalizar solo MUEVE (Move-Item): nunca borra ni sobrescribe. Si el destino ya tiene un archivo con ese
nombre, lo deja donde esta y lo anota como CONFLICTO. Todo queda en normalizacion_log.csv. Las carpetas que
queden vacias no se borran.
Salida: pipeline_integridad.csv (una fila por .tif original) + resumen en pantalla.
Con -Tabla: RawData_actualizado.csv (misma tabla y mismo orden de filas; el original NO se toca).
  Reglas por columna (limpieza / Stardist / ROI analysis), segun la verificacion de ese paso:
    valido y completo -> "ok"        con error o incoherente -> "revisar"
    no existe todavia -> la celda queda como estaba (no se borran tus notas)
  Emparejamiento de filas: por nombre de imagen (sin .tif ni sufijo _Stack/_ROI) Y por mutante: si la
  tabla tiene la columna "mut", solo se consideran las carpetas cuyo primer nivel bajo -Raiz se llama
  igual (ej. ...\\1265\\...). Si el nombre no coincide pero el "abf file" de la fila coincide con el .abf
  que hay junto al .tif, se empareja por el .abf y se avisa (para corregir el nombre).
  Completa "mut" y "abf file" solo cuando estan vacios. Los .tif del disco que no estan en la tabla se
  AGREGAN al final como filas nuevas (mut = carpeta de primer nivel, abf = el .abf de la carpeta si hay
  uno solo); placa/corte/region/obs quedan vacios. Para no agregar filas usa -NoAgregar.
  Despues copia las columnas y pegalas en la hoja.
Si Windows bloquea el script:
  powershell -ExecutionPolicy Bypass -File .\verificar_pipeline.ps1 -Raiz "D:\..."
#>
param(
    [Parameter(Mandatory = $true)][string]$Raiz,
    [string]$Salida = "pipeline_integridad.csv",
    [string]$CarpetaLimpieza = "limpieza_ij",
    [string]$CarpetaRois = "rois_ij",
    [string[]]$ExcluirCarpetas = @("ROIs", "limpieza_ij", "procesado", "batch_analysis", "rois_ij"),
    [switch]$SinROIs,
    [switch]$Normalizar,
    [switch]$Aplicar,
    [string]$Tabla = "",
    [string]$SalidaTabla = "RawData_actualizado.csv",
    [string]$Marca = "ok",
    [string]$MarcaRevisar = "revisar",
    [switch]$NoAgregar,
    [string]$ColMut = "mut",
    [string]$ColAbf = "abf file",
    [string]$ColImagen = "image File",
    [string]$ColLimpieza = "limpieza",
    [string]$ColStardist = "Stardist",
    [string]$ColRoi = "ROI analysis"
)

$ErrorActionPreference = "Stop"
if ($Aplicar -and -not $Normalizar) {
    Write-Host "-Aplicar solo se usa junto con -Normalizar." -ForegroundColor Red
    exit 1
}
Add-Type -AssemblyName System.IO.Compression.FileSystem

if (-not (Test-Path -LiteralPath $Raiz)) {
    Write-Host "No existe la carpeta raiz: $Raiz" -ForegroundColor Red
    exit 1
}

$RaizAbs = (Resolve-Path -LiteralPath $Raiz).ProviderPath.TrimEnd('\', '/')

# ---------- utilidades ----------
function Get-TiffFrames([string]$ruta) {
    # Lee los primeros 64 KB y busca "images=N" de la descripcion ImageJ. Devuelve -1 si no hay.
    try {
        $fs = [IO.File]::Open($ruta, 'Open', 'Read', 'ReadWrite')
        try {
            $n = [int][Math]::Min(65536, $fs.Length)
            $buf = New-Object byte[] $n
            [void]$fs.Read($buf, 0, $n)
        } finally { $fs.Close() }
        $txt = [Text.Encoding]::GetEncoding(28591).GetString($buf)
        $m = [regex]::Match($txt, 'images=(\d+)')
        if ($m.Success) { return [int]$m.Groups[1].Value }
        return -1
    } catch { return -2 }
}

function Get-LineCount([string]$ruta) {
    $c = 0
    $sr = New-Object IO.StreamReader($ruta)
    try { while ($null -ne $sr.ReadLine()) { $c++ } } finally { $sr.Close() }
    return $c
}

function Get-RoiCountZip([string]$ruta) {
    # -1 = no se pudo abrir como zip
    try {
        $z = [IO.Compression.ZipFile]::OpenRead($ruta)
        try {
            $n = 0
            foreach ($e in $z.Entries) { if ($e.Name -match '\.roi$') { $n++ } }
            return $n
        } finally { $z.Dispose() }
    } catch { return -1 }
}

function Get-AreaCols([string]$ruta) {
    $sr = New-Object IO.StreamReader($ruta)
    try { $h = $sr.ReadLine() } finally { $sr.Close() }
    if ($null -eq $h) { return -1 }
    $n = 0
    foreach ($c in ($h -split ',')) { if ($c.Trim('"').Trim() -match '^Area\d+$') { $n++ } }
    return $n
}

function Get-Clave([string]$nombre) {
    $t = $nombre.Trim().ToLower()
    $t = $t -replace '\.(tif|tiff)$', ''
    $t = $t -replace '_(stack|roi)$', ''
    return $t
}

# ---------- 1. localizar los .tif originales ----------
Write-Host "Buscando archivos en $Raiz ..."
$todosArch = Get-ChildItem -LiteralPath $Raiz -Recurse -File -ErrorAction SilentlyContinue
$todos = $todosArch | Where-Object { $_.Extension -match '^\.tiff?$' }

# Indice de salidas del pipeline por nombre de archivo (en cualquier subcarpeta)
$idx = @{}
foreach ($a in $todosArch) {
    if ($a.Name -match '_clean_quant\.tif+$|_clean_view\.tif+$|_drift_shifts_ij\.csv$|_RoiSet\.zip$|_ROI_Results\.csv$|_rois_table_ij\.csv$|_rois_labels\.tif+$') {
        $kk = $a.Name.ToLower()
        if (-not $idx.ContainsKey($kk)) { $idx[$kk] = New-Object System.Collections.ArrayList }
        [void]$idx[$kk].Add($a.FullName)
    }
}

# Archivos .abf por carpeta
$abfPorDir = @{}
foreach ($a in $todosArch) {
    if ($a.Extension -ieq ".abf") {
        $kd = $a.DirectoryName.ToLower()
        if (-not $abfPorDir.ContainsKey($kd)) { $abfPorDir[$kd] = New-Object System.Collections.ArrayList }
        [void]$abfPorDir[$kd].Add($a.Name)
    }
}

# Busca una salida por nombre: primero bajo la carpeta del original (cualquier subcarpeta), luego en la
# carpeta superior (mismo nivel). Devuelve $fallback (ruta inexistente) si no la encuentra.
function Find-Salida([string]$nombre, [string]$dirRaw, [string]$fallback) {
    $kk = $nombre.ToLower()
    if (-not $idx.ContainsKey($kk)) { return $fallback }
    $base = $dirRaw.TrimEnd('\', '/') + '\'
    $padre = Split-Path $dirRaw -Parent
    if ($null -eq $padre -or $padre -eq "") { $padre = $dirRaw }
    $padre = $padre.TrimEnd('\', '/') + '\'
    $mejor = $null; $mejorScore = 0; $mejorLen = 0
    foreach ($p in $idx[$kk]) {
        $d = (Split-Path $p -Parent).TrimEnd('\', '/') + '\'
        $score = 0
        if ($d.StartsWith($base, [StringComparison]::OrdinalIgnoreCase)) { $score = 2 }
        elseif ($d.StartsWith($padre, [StringComparison]::OrdinalIgnoreCase)) { $score = 1 }
        if ($score -eq 0) { continue }
        if ($score -gt $mejorScore -or ($score -eq $mejorScore -and $p.Length -lt $mejorLen)) {
            $mejor = $p; $mejorScore = $score; $mejorLen = $p.Length
        }
    }
    if ($null -eq $mejor) { return $fallback }
    return $mejor
}

$originales = @()
foreach ($f in $todos) {
    if ($f.Name -match '_clean_(quant|view)\.tif+$|_rois_labels\.tif+$') { continue }
    $rel = $f.DirectoryName.Substring([Math]::Min($RaizAbs.Length, $f.DirectoryName.Length))
    $segs = $rel -split '[\\/]' | Where-Object { $_ -ne "" }
    $excl = $false
    foreach ($s in $segs) {
        foreach ($x in $ExcluirCarpetas) { if ($s -ieq $x) { $excl = $true } }
    }
    if (-not $excl) { $originales += $f }
}
Write-Host ("TIFF originales encontrados: {0}" -f $originales.Count)
if ($originales.Count -eq 0) {
    Write-Host "No se encontro ninguno. Revisa -Raiz o -ExcluirCarpetas." -ForegroundColor Yellow
    exit 1
}

$logNorm = New-Object System.Collections.ArrayList
function Normalizar-Archivo([string]$ruta, [string]$dirEsp, [string]$stem) {
    # Mueve $ruta a $dirEsp si esta en otra carpeta. Devuelve la ruta resultante.
    if (-not (Test-Path -LiteralPath $ruta)) { return $ruta }
    $dAct = (Split-Path $ruta -Parent).TrimEnd('\', '/')
    if ($dAct -ieq $dirEsp.TrimEnd('\', '/')) { return $ruta }
    $leaf = Split-Path $ruta -Leaf
    $dest = Join-Path $dirEsp $leaf
    $reg = [PSCustomObject]@{ Stack = $stem; Archivo = $leaf; Desde = $dAct; Hacia = $dirEsp; Accion = "" }
    if (Test-Path -LiteralPath $dest) {
        $reg.Accion = "CONFLICTO (ya existe en destino; no se movio)"
        [void]$logNorm.Add($reg)
        return $ruta
    }
    if ($Aplicar) {
        try {
            if (-not (Test-Path -LiteralPath $dirEsp)) { [void](New-Item -ItemType Directory -Path $dirEsp -Force) }
            Move-Item -LiteralPath $ruta -Destination $dest
            $reg.Accion = "MOVIDO"
            [void]$logNorm.Add($reg)
            return $dest
        } catch {
            $reg.Accion = "ERROR: " + $_.Exception.Message
            [void]$logNorm.Add($reg)
            return $ruta
        }
    }
    $reg.Accion = "PLAN (no aplicado)"
    [void]$logNorm.Add($reg)
    return $ruta
}

# ---------- 2. verificar cada uno ----------
$res = @()
$i = 0
foreach ($f in $originales) {
    $i++
    if (($i % 10) -eq 0) { Write-Host ("  {0}/{1}" -f $i, $originales.Count) }
    $dir = $f.DirectoryName
    $stem = [IO.Path]::GetFileNameWithoutExtension($f.Name)
    $dl = Join-Path $dir $CarpetaLimpieza
    $dr = Join-Path $dir $CarpetaRois

    $pQuant = Find-Salida ($stem + "_clean_quant.tif") $dir (Join-Path $dl ($stem + "_clean_quant.tif"))
    $pView  = Find-Salida ($stem + "_clean_view.tif") $dir (Join-Path $dl ($stem + "_clean_view.tif"))
    $pDrift = Find-Salida ($stem + "_drift_shifts_ij.csv") $dir (Join-Path $dl ($stem + "_drift_shifts_ij.csv"))
    $pZip   = Find-Salida ($stem + "_RoiSet.zip") $dir (Join-Path $dr ($stem + "_RoiSet.zip"))
    $pCsv   = Find-Salida ($stem + "_ROI_Results.csv") $dir (Join-Path $dr ($stem + "_ROI_Results.csv"))

    $pTabla  = Find-Salida ($stem + "_rois_table_ij.csv") $dir (Join-Path $dr ($stem + "_rois_table_ij.csv"))
    $pLabels = Find-Salida ($stem + "_rois_labels.tif") $dir (Join-Path $dl ($stem + "_rois_labels.tif"))
    if ($Normalizar) {
        $pQuant  = Normalizar-Archivo $pQuant $dl $stem
        $pView   = Normalizar-Archivo $pView $dl $stem
        $pDrift  = Normalizar-Archivo $pDrift $dl $stem
        $pLabels = Normalizar-Archivo $pLabels $dl $stem
        $pZip    = Normalizar-Archivo $pZip $dr $stem
        $pCsv    = Normalizar-Archivo $pCsv $dr $stem
        $pTabla  = Normalizar-Archivo $pTabla $dr $stem
    }

    $errL = New-Object System.Collections.ArrayList
    $errS = New-Object System.Collections.ArrayList
    $errR = New-Object System.Collections.ArrayList
    $avi = New-Object System.Collections.ArrayList
    $estado = "OK"

    $hayQuant = Test-Path -LiteralPath $pQuant
    $hayZip = Test-Path -LiteralPath $pZip
    $hayCsv = Test-Path -LiteralPath $pCsv

    $fr = Get-TiffFrames $f.FullName
    $frQ = -1; $nRoiZip = -1; $nRoiCsv = -1; $filasCsv = -1; $filasDrift = -1

    # --- paso 1
    if (-not $hayQuant) {
        $estado = "FALTA_LIMPIEZA"
    } else {
        $q = Get-Item -LiteralPath $pQuant
        if ($q.Length -eq 0) { [void]$errL.Add("clean_quant vacio (0 bytes)") }
        else {
            $frQ = Get-TiffFrames $pQuant
            if ($fr -ge 0 -and $frQ -ge 0 -and $fr -ne $frQ) { [void]$errL.Add("frames original=$fr, clean_quant=$frQ") }
            if ($q.LastWriteTime -lt $f.LastWriteTime) { [void]$avi.Add("clean_quant mas antiguo que el original") }
        }
        if (-not (Test-Path -LiteralPath $pView)) { [void]$avi.Add("falta clean_view") }
        elseif ((Get-Item -LiteralPath $pView).Length -eq 0) { [void]$errL.Add("clean_view vacio") }
        if (-not (Test-Path -LiteralPath $pDrift)) { [void]$avi.Add("falta CSV de deriva") }
        else {
            $filasDrift = (Get-LineCount $pDrift) - 1
            if ($filasDrift -le 0) { [void]$errL.Add("CSV de deriva vacio") }
            elseif ($frQ -ge 0 -and $filasDrift -ne $frQ) { [void]$errL.Add("CSV de deriva: $filasDrift filas, clean_quant: $frQ frames") }
        }
    }

    # --- pasos 2 y 3
    if (-not $SinROIs -and $hayQuant) {
        if (-not $hayZip) {
            if ($estado -eq "OK") { $estado = "FALTA_ROIS" }
        } else {
            $z = Get-Item -LiteralPath $pZip
            if ($z.Length -eq 0) { [void]$errS.Add("RoiSet vacio (0 bytes)") }
            else {
                $nRoiZip = Get-RoiCountZip $pZip
                if ($nRoiZip -lt 0) { [void]$errS.Add("RoiSet.zip corrupto (no se abre)") }
                elseif ($nRoiZip -eq 0) { [void]$errS.Add("RoiSet sin ROIs") }
            }
            if ((Test-Path -LiteralPath $pQuant) -and ($z.LastWriteTime -lt (Get-Item -LiteralPath $pQuant).LastWriteTime)) {
                [void]$avi.Add("RoiSet mas antiguo que clean_quant")
            }
        }
        if (-not $hayCsv) {
            if ($hayZip -and $estado -eq "OK") { $estado = "FALTA_MEDICION" }
        } else {
            $c = Get-Item -LiteralPath $pCsv
            if ($c.Length -eq 0) { [void]$errR.Add("ROI_Results vacio (0 bytes)") }
            else {
                $nRoiCsv = Get-AreaCols $pCsv
                $filasCsv = (Get-LineCount $pCsv) - 1
                if ($nRoiCsv -le 0) { [void]$errR.Add("ROI_Results sin columnas Area{k}") }
                if ($nRoiZip -gt 0 -and $nRoiCsv -gt 0 -and $nRoiZip -ne $nRoiCsv) {
                    [void]$errR.Add("ROI: RoiSet=$nRoiZip, ROI_Results=$nRoiCsv")
                }
                if ($frQ -ge 0 -and $filasCsv -ne $frQ) { [void]$errR.Add("ROI_Results: $filasCsv filas, clean_quant: $frQ frames") }
                if ($hayZip -and ($c.LastWriteTime -lt (Get-Item -LiteralPath $pZip).LastWriteTime)) {
                    [void]$avi.Add("ROI_Results mas antiguo que el RoiSet (medir de nuevo?)")
                }
            }
        }
    }

    # ubicacion no estandar (no es error, pero conviene saberlo)
    $esperadas = @(
        @($pQuant, $dl), @($pView, $dl), @($pDrift, $dl), @($pZip, $dr), @($pCsv, $dr)
    )
    foreach ($par in $esperadas) {
        if (Test-Path -LiteralPath $par[0]) {
            $dd = (Split-Path $par[0] -Parent).TrimEnd('\', '/')
            if ($dd -ine $par[1].TrimEnd('\', '/')) {
                [void]$avi.Add(("{0} fuera de su carpeta esperada: {1}" -f (Split-Path $par[0] -Leaf), $dd))
            }
        }
    }

    $err = @($errL) + @($errS) + @($errR)
    if ($err.Count -gt 0) { $estado = "ERROR" }

    # estado por paso: ok / error / falta
    $stL = "falta"; $stS = "falta"; $stR = "falta"
    if ($hayQuant) { if ($errL.Count -eq 0) { $stL = "ok" } else { $stL = "error" } }
    if ($SinROIs) { $stS = "n/a"; $stR = "n/a" }
    else {
        if ($hayQuant -and $hayZip) { if ($errS.Count -eq 0) { $stS = "ok" } else { $stS = "error" } }
        if ($hayQuant -and $hayCsv) { if ($errR.Count -eq 0) { $stR = "ok" } else { $stR = "error" } }
    }

    $relDir = $dir.Substring([Math]::Min($RaizAbs.Length, $dir.Length)).Trim('\', '/')
    $mutDe = ""
    if ($relDir -ne "") { $mutDe = ($relDir -split '[\\/]')[0] }
    $abfs = @()
    if ($abfPorDir.ContainsKey($dir.ToLower())) { $abfs = @($abfPorDir[$dir.ToLower()]) }
    $abfKeys = @($abfs | ForEach-Object { [IO.Path]::GetFileNameWithoutExtension($_).ToLower() })

    $res += [PSCustomObject]@{
        Estado        = $estado
        Stack         = $stem
        Clave         = (Get-Clave $stem)
        Mut           = $mutDe
        Abf           = ($abfs -join ";")
        AbfKeys       = $abfKeys
        Carpeta       = $dir
        MB            = [Math]::Round($f.Length / 1MB, 1)
        Frames_orig   = $fr
        Frames_clean  = $frQ
        Filas_deriva  = $filasDrift
        ROIs_RoiSet   = $nRoiZip
        ROIs_Results  = $nRoiCsv
        Filas_Results = $filasCsv
        Limpieza      = $stL
        Stardist      = $stS
        ROI_analysis  = $stR
        Detalle       = ($err -join " | ")
        Avisos        = ($avi -join " | ")
    }
}

# ---------- 3. salida ----------
$res | Sort-Object Estado, Carpeta, Stack | Select-Object -Property * -ExcludeProperty Clave, AbfKeys |
    Export-Csv -Path $Salida -NoTypeInformation -Encoding UTF8 -Delimiter ","

Write-Host ""
Write-Host "=== RESUMEN ===" -ForegroundColor Cyan
$res | Group-Object Estado | Sort-Object Name | ForEach-Object {
    Write-Host ("  {0,-16} {1,4}" -f $_.Name, $_.Count)
}
$malos = $res | Where-Object { $_.Estado -ne "OK" }
if (@($malos).Count -gt 0) {
    Write-Host ""
    Write-Host "Stacks que requieren atencion:" -ForegroundColor Yellow
    foreach ($m in ($malos | Sort-Object Estado, Stack)) {
        $t = "  [{0}] {1}" -f $m.Estado, $m.Stack
        if ($m.Detalle -ne "") { $t = $t + "  -> " + $m.Detalle }
        Write-Host $t
    }
}
$conAvisos = $res | Where-Object { $_.Avisos -ne "" }
if (@($conAvisos).Count -gt 0) {
    Write-Host ""
    Write-Host ("Stacks con avisos (no son errores): {0}  (ver columna Avisos)" -f @($conAvisos).Count) -ForegroundColor DarkYellow
}
Write-Host ""
Write-Host "Detalle completo en: $Salida"
Write-Host "Nota: nombres de stack repetidos en distintas carpetas se verifican por separado (cada uno en su carpeta)."


if ($Normalizar) {
    Write-Host ""
    Write-Host "=== NORMALIZACION DE UBICACIONES ===" -ForegroundColor Cyan
    if ($logNorm.Count -eq 0) {
        Write-Host "Todos los archivos ya estaban en su carpeta estandar."
    } else {
        $logNorm | Export-Csv -Path "normalizacion_log.csv" -NoTypeInformation -Encoding UTF8 -Delimiter ","
        $logNorm | Group-Object { ($_.Accion -split ' ')[0] } | ForEach-Object { Write-Host ("  {0,-12} {1,4}" -f $_.Name, $_.Count) }
        if (-not $Aplicar) {
            Write-Host "Era solo una simulacion: no se movio nada. Revisa normalizacion_log.csv y, si esta bien," -ForegroundColor Yellow
            Write-Host "repite el comando agregando -Aplicar." -ForegroundColor Yellow
        }
        Write-Host "Detalle en: normalizacion_log.csv"
    }
}

# ---------- 4. actualizar la tabla Raw Data (opcional) ----------
if ($Tabla -ne "") {
    if (-not (Test-Path -LiteralPath $Tabla)) {
        Write-Host "No existe la tabla: $Tabla" -ForegroundColor Red
        exit 1
    }
    # Delimiter explicito: en Windows en espanol el predeterminado seria ";"
    $filas = @(Import-Csv -LiteralPath $Tabla -Delimiter ",")
    if ($filas.Count -eq 0) { Write-Host "La tabla esta vacia." -ForegroundColor Red; exit 1 }
    $cols = $filas[0].PSObject.Properties.Name
    $necesarias = @($ColImagen, $ColLimpieza)
    if (-not $SinROIs) { $necesarias += @($ColStardist, $ColRoi) }
    foreach ($c in $necesarias) {
        if ($cols -notcontains $c) {
            Write-Host "Falta la columna '$c' en la tabla. Columnas encontradas: $($cols -join ' | ')" -ForegroundColor Red
            exit 1
        }
    }
    $tieneMut = ($cols -contains $ColMut)
    $tieneAbf = ($cols -contains $ColAbf)
    if (-not $tieneMut) { Write-Host "Aviso: la tabla no tiene la columna '$ColMut'; se empareja solo por nombre." -ForegroundColor Yellow }

    $porClave = @{}
    foreach ($r in $res) {
        if (-not $porClave.ContainsKey($r.Clave)) { $porClave[$r.Clave] = New-Object System.Collections.ArrayList }
        [void]$porClave[$r.Clave].Add($r)
    }

    function Mut-Coincide([string]$a, [string]$b) {
        if ($a.Trim() -eq "" -or $b.Trim() -eq "") { return $true }
        return ($a.Trim() -ieq $b.Trim())
    }
    function Set-Celda($fila, [string]$col, [string[]]$estados) {
        # devuelve 1 si cambio la celda
        $nuevo = $null
        if ($estados -contains "error") { $nuevo = $MarcaRevisar }
        elseif ($estados.Count -gt 0 -and (($estados | Where-Object { $_ -ne "ok" }).Count -eq 0)) { $nuevo = $Marca }
        if ($null -ne $nuevo -and $fila.$col -ne $nuevo) { $fila.$col = $nuevo; return 1 }
        return 0
    }
    function Valor-Estado([string]$st) {
        if ($st -eq "ok") { return $Marca }
        if ($st -eq "error") { return $MarcaRevisar }
        return ""
    }

    $cambios = 0
    $usados = @{}
    $sinDisco = @()
    $porAbf = @()
    foreach ($r in $filas) {
        $nombre = [string]$r.$ColImagen
        if ($nombre.Trim() -eq "") { continue }
        $k = Get-Clave $nombre
        $filaMut = ""
        if ($tieneMut) { $filaMut = [string]$r.$ColMut }
        $m = @()
        if ($porClave.ContainsKey($k)) { $m = @($porClave[$k] | Where-Object { Mut-Coincide $filaMut $_.Mut }) }
        $via = "nombre"
        if ($m.Count -eq 0 -and $tieneAbf) {
            $ak = [IO.Path]::GetFileNameWithoutExtension(([string]$r.$ColAbf).Trim()).ToLower()
            if ($ak -ne "") {
                $m = @($res | Where-Object { ($_.AbfKeys -contains $ak) -and (Mut-Coincide $filaMut $_.Mut) })
                if ($m.Count -gt 0) { $via = "abf" }
            }
        }
        if ($m.Count -eq 0) { $sinDisco += $nombre; continue }
        foreach ($x in $m) { $usados[$x.Carpeta + "|" + $x.Stack] = $true }
        if ($via -eq "abf") { $porAbf += ("tabla: {0}  ->  disco: {1}  ({2})" -f $nombre, $m[0].Stack, $m[0].Carpeta) }

        $cambios += Set-Celda $r $ColLimpieza @($m | ForEach-Object { $_.Limpieza })
        if (-not $SinROIs) {
            $cambios += Set-Celda $r $ColStardist @($m | ForEach-Object { $_.Stardist })
            $cambios += Set-Celda $r $ColRoi @($m | ForEach-Object { $_.ROI_analysis })
        }
        # completar mut y abf solo si estan vacios y la coincidencia es unica
        if ($m.Count -eq 1) {
            if ($tieneMut -and ([string]$r.$ColMut).Trim() -eq "" -and $m[0].Mut -ne "") { $r.$ColMut = $m[0].Mut; $cambios++ }
            if ($tieneAbf -and ([string]$r.$ColAbf).Trim() -eq "" -and $m[0].Abf -ne "" -and $m[0].Abf -notmatch ';') { $r.$ColAbf = $m[0].Abf; $cambios++ }
        }
    }

    # TIFF en disco que no estan en la tabla -> filas nuevas al final
    $nuevos = @($res | Where-Object { -not $usados.ContainsKey($_.Carpeta + "|" + $_.Stack) } | Sort-Object Mut, Carpeta, Stack)
    $agregadas = @()
    if ($nuevos.Count -gt 0 -and -not $NoAgregar) {
        foreach ($n in $nuevos) {
            $o = [ordered]@{}
            foreach ($c in $cols) { $o[$c] = "" }
            $o[$ColImagen] = $n.Stack
            if ($tieneMut) { $o[$ColMut] = $n.Mut }
            if ($tieneAbf -and $n.Abf -ne "" -and $n.Abf -notmatch ';') { $o[$ColAbf] = $n.Abf }
            $o[$ColLimpieza] = Valor-Estado $n.Limpieza
            if (-not $SinROIs) {
                $o[$ColStardist] = Valor-Estado $n.Stardist
                $o[$ColRoi] = Valor-Estado $n.ROI_analysis
            }
            $agregadas += [PSCustomObject]$o
        }
        $filas = @($filas) + $agregadas
    }
    $filas | Export-Csv -Path $SalidaTabla -NoTypeInformation -Encoding UTF8 -Delimiter ","

    Write-Host ""
    Write-Host "=== TABLA RAW DATA ===" -ForegroundColor Cyan
    Write-Host ("Celdas modificadas: {0}  |  filas nuevas agregadas: {1}  -> {2}" -f $cambios, $agregadas.Count, $SalidaTabla) -ForegroundColor Green
    if ($agregadas.Count -gt 0) {
        Write-Host "Filas nuevas (al final de la tabla; faltan placa / corte / region / obs):" -ForegroundColor Yellow
        $agregadas | ForEach-Object { Write-Host ("  {0}  [{1}]" -f $_.$ColImagen, $(if ($tieneMut) { $_.$ColMut } else { "" })) }
    }
    if ($nuevos.Count -gt 0 -and $NoAgregar) {
        Write-Host "TIFF en disco que NO estan en la tabla (no se agregaron por -NoAgregar):" -ForegroundColor Yellow
        $nuevos | ForEach-Object { Write-Host "  $($_.Stack)   ($($_.Carpeta))" }
    }
    if ($porAbf.Count -gt 0) {
        Write-Host ""
        Write-Host "Emparejadas por archivo .abf porque el NOMBRE de imagen no coincide (revisar cual es el correcto):" -ForegroundColor Yellow
        $porAbf | ForEach-Object { Write-Host "  $_" }
    }
    if ($sinDisco.Count -gt 0) {
        Write-Host ""
        Write-Host "Filas de la tabla sin ningun .tif en disco (celdas sin tocar):" -ForegroundColor Yellow
        $sinDisco | Sort-Object -Unique | ForEach-Object { Write-Host "  $_" }
    }
}

if (@($malos).Count -eq 0) { exit 0 } else { exit 2 }

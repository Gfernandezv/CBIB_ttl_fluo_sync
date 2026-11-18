"""
diagnostico_abf.py - Revisa los .abf de TODAS las muestras del pipeline y los compara con ROI_Results.

Por muestra reporta: frecuencia de muestreo, rango de temperatura y TTL, nº de pulsos TTL detectados,
nº de frames de ROI_Results, intervalo entre pulsos (mediana y variacion) y variacion de la amplitud
de los pulsos. Sirve para ver si la amplitud irregular del TTL y/o el desajuste pulsos-frames se
repite en todas las grabaciones o es de una sola.

Uso (desde la carpeta scripts/ o con la ruta completa):
    python diagnostico_abf.py "/Volumes/TOSHIBA EXT/Germán Fernandez/TRPM3/Fluo4"
    python diagnostico_abf.py RAIZ --threshold 0.1 --ttl-channel 1 --temp-channel 0 --salida diag.csv

Solo lee; no modifica nada. Escribe un CSV (por defecto diagnostico_abf.csv en la carpeta actual).
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyabf


def pulses(ttl, t, thr):
    b = ttl > thr
    d = np.diff(b.astype(int))
    on = np.where(d == 1)[0] + 1
    off = np.where(d == -1)[0] + 1
    if b[0]:
        on = np.r_[0, on]
    if b[-1]:
        off = np.r_[off, len(b) - 1]
    n = min(len(on), len(off))
    return on[:n], off[:n]


def diagnose(sample_dir, roi_csv, abf_path, thr, ttl_ch, temp_ch):
    abf = pyabf.ABF(str(abf_path))
    abf.setSweep(0, channel=temp_ch); t = abf.sweepX.copy(); temp = abf.sweepY.copy()
    abf.setSweep(0, channel=ttl_ch); ttl = abf.sweepY.copy()
    on, off = pulses(ttl, t, thr)
    n_frames = len(pd.read_csv(roi_csv))
    row = {
        "abf": abf_path.name,
        "n_abf_en_carpeta": len(list(Path(sample_dir).glob("*.abf"))),
        "canales": "|".join(abf.adcNames),
        "fs_hz": abf.dataRate,
        "dur_s": round(float(abf.sweepLengthSec), 1),
        "temp_min": round(float(np.nanmin(temp)), 1),
        "temp_max": round(float(np.nanmax(temp)), 1),
        "ttl_min": round(float(np.nanmin(ttl)), 3),
        "ttl_max": round(float(np.nanmax(ttl)), 3),
        "n_pulsos": len(on),
        "n_frames": n_frames,
        "pulsos_menos_frames": len(on) - n_frames,
    }
    if len(on) > 2:
        d = np.diff(t[on])
        peaks = np.array([ttl[a:max(b, a + 1)].max() for a, b in zip(on, off)])
        width = (t[off] - t[on])
        row.update({
            "intervalo_med_s": round(float(np.median(d)), 4),
            "intervalo_cv": round(float(np.std(d) / np.mean(d)), 3),
            "intervalo_max_s": round(float(d.max()), 3),
            "ancho_pulso_med_ms": round(float(np.median(width) * 1000), 2),
            "muestras_por_pulso": round(float(np.median(width) * abf.dataRate), 1),
            "amp_pulso_med": round(float(np.median(peaks)), 3),
            "amp_pulso_cv": round(float(np.std(peaks) / np.mean(peaks)), 3),
            "amp_pulso_min": round(float(peaks.min()), 3),
        })
    flags = []
    if len(on) == 0:
        flags.append("SIN_PULSOS")
    elif len(on) < n_frames:
        flags.append("FALTAN_PULSOS")
    elif len(on) > n_frames * 1.05:
        flags.append("SOBRAN_PULSOS")
    if row.get("amp_pulso_cv", 0) > 0.2:
        flags.append("AMPLITUD_IRREGULAR")
    if row.get("muestras_por_pulso", 99) < 3:
        flags.append("PULSO_CORTO_vs_fs")
    if row["n_abf_en_carpeta"] > 1:
        flags.append("VARIOS_ABF")
    row["flags"] = " ".join(flags) or "ok"
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raiz")
    ap.add_argument("--threshold", type=float, default=0.1)
    ap.add_argument("--ttl-channel", type=int, default=1)
    ap.add_argument("--temp-channel", type=int, default=0)
    ap.add_argument("--salida", default="diagnostico_abf.csv")
    a = ap.parse_args()
    raiz = Path(a.raiz)
    rows = []
    for roi_csv in sorted(raiz.rglob("*_ROI_Results.csv")):
        d = roi_csv.parent
        while d != raiz and d.name.lower() in {"rois_ij", "limpieza_ij"}:
            d = d.parent
        stem = roi_csv.name[: -len("_ROI_Results.csv")]
        base = {"folder": d.relative_to(raiz).as_posix(), "stem": stem}
        abfs = sorted(d.glob("*.abf"))
        if not abfs:
            rows.append({**base, "flags": "SIN_ABF"}); continue
        try:
            rows.append({**base, **diagnose(d, roi_csv, abfs[0], a.threshold, a.ttl_channel, a.temp_channel)})
        except Exception as e:
            rows.append({**base, "abf": abfs[0].name, "flags": f"ERROR {type(e).__name__}: {e}"})
    df = pd.DataFrame(rows)
    df.to_csv(a.salida, index=False)
    cols = [c for c in ["folder", "abf", "fs_hz", "ttl_max", "n_pulsos", "n_frames", "pulsos_menos_frames",
                        "intervalo_cv", "amp_pulso_cv", "muestras_por_pulso", "flags"] if c in df.columns]
    pd.set_option("display.width", 250, "display.max_colwidth", 60, "display.max_rows", 500)
    print(df[cols].to_string(index=False))
    print(f"\n{len(df)} muestras | detalle: {a.salida}")
    print(df["flags"].str.split().explode().value_counts().to_string())


if __name__ == "__main__":
    main()

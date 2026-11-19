"""Registro (manifiesto) de datos: recorre <raiz>/<genotipo>/**/ y deja UNA tabla con cada registro, su umbral TTL y
si entra al analisis o se excluye (con el motivo). El nombre de la carpeta de primer nivel es el genotipo."""
import os
import re
import glob
import numpy as np
import pandas as pd
from . import signal as sg

COLUMNAS = ["genotipo", "registro", "carpeta", "abf", "roi_csv", "deriva_csv", "frames", "pulsos", "umbral_ttl", "incluido", "motivo"]


def _numero_registro(rc, carpeta):
    """Numero del registro a partir del nombre del archivo de ROI. Acepta <fecha>_<n>_Stack_ROI_Results.csv y
    <fecha>_<n>_ROI_ROI_Results.csv; si no hay numero, usa el nombre de la carpeta de la muestra."""
    m = re.search(r"_(\d+)(?:_|$)", os.path.basename(rc))
    return int(m.group(1)) if m else os.path.basename(carpeta)


def escanear(raiz, cfg):
    filas = []
    ttl = cfg["ttl"]
    for geno in cfg["genotipos"]:
        patron = os.path.join(raiz, geno, "**", "rois_ij", "*_ROI_Results.csv")
        for rc in sorted(glob.glob(patron, recursive=True)):
            carpeta = os.path.dirname(os.path.dirname(rc))
            reg = _numero_registro(rc, carpeta)
            abfs = glob.glob(os.path.join(carpeta, "*.abf"))
            deriva = glob.glob(os.path.join(carpeta, "limpieza_ij", "*drift_shifts_ij.csv"))
            fila = dict(genotipo=geno, registro=reg, carpeta=os.path.relpath(carpeta, raiz), abf="", roi_csv=os.path.relpath(rc, raiz),
                        deriva_csv="", frames=np.nan, pulsos=np.nan, umbral_ttl=np.nan, incluido=False, motivo="")
            if len(abfs) != 1 or not deriva:
                fila["motivo"] = "falta el abf (o hay mas de uno) o la tabla de deriva"; filas.append(fila); continue
            fila["abf"] = os.path.relpath(abfs[0], raiz); fila["deriva_csv"] = os.path.relpath(deriva[0], raiz)
            n = int(pd.read_csv(rc, usecols=[0]).shape[0]); fila["frames"] = n
            _, _, v = sg.leer_abf(abfs[0], ttl["canal_temperatura"], ttl["canal_ttl"])
            fijo = (ttl.get("umbral_fijo") or {}).get(geno)
            umbral = fijo if fijo is not None else (sg.buscar_umbral(v, n) if ttl["umbral"] == "auto" else float(ttl["umbral"]))
            if umbral is None:
                umbral = float((np.median(v) + v.max()) / 2)
            fila["umbral_ttl"] = round(float(umbral), 4); fila["pulsos"] = int(len(sg.pulsos(v, umbral)[0]))
            if abs(fila["pulsos"] - n) > ttl["tolerancia_pulsos"]:
                fila["motivo"] = f"pulsos TTL ({fila['pulsos']}) no coinciden con frames ({n})"
            else:
                fila["incluido"] = True
            filas.append(fila)
    return pd.DataFrame(filas, columns=COLUMNAS)


def detalle_pulsos(registro, raiz, cfg):
    """Tabla de diagnostico de pulsos TTL por registro (con abf), para entender las exclusiones."""
    ttl = cfg["ttl"]; filas = []
    for _, f in registro[registro.abf != ""].iterrows():
        t, _, v = sg.leer_abf(os.path.join(raiz, f["abf"]), ttl["canal_temperatura"], ttl["canal_ttl"])
        d = sg.diagnostico_pulsos(t, v, f["umbral_ttl"], int(f["frames"]))
        filas.append(dict(genotipo=f["genotipo"], registro=f["registro"], umbral_ttl=f["umbral_ttl"], incluido=f["incluido"], **d))
    return pd.DataFrame(filas)

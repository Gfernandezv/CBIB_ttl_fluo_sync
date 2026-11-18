"""Registro (manifiesto) de datos: recorre <raiz>/<genotipo>/**/ y deja UNA tabla con cada registro, su umbral TTL y
si entra al analisis o se excluye (con el motivo). El nombre de la carpeta de primer nivel es el genotipo."""
import os
import glob
import numpy as np
import pandas as pd
from . import signal as sg

COLUMNAS = ["genotipo", "registro", "carpeta", "abf", "roi_csv", "deriva_csv", "frames", "pulsos", "umbral_ttl", "incluido", "motivo"]


def escanear(raiz, cfg):
    filas = []
    ttl = cfg["ttl"]
    for geno in cfg["genotipos"]:
        patron = os.path.join(raiz, geno, "**", "rois_ij", "*_Stack_ROI_Results.csv")
        for rc in sorted(glob.glob(patron, recursive=True)):
            carpeta = os.path.dirname(os.path.dirname(rc))
            reg = int(os.path.basename(rc).split("_")[1])
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

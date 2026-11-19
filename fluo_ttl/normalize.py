"""ΔF/F0 corregido por la deriva comun de cada registro.

F0 (por ROI)  = media de los `frames_f0` primeros frames, antes de calentar.
c(t)          = deriva comun: mediana entre ROIs de F/F0 SOLO en frames basales (T <= T0 + margen, sin transicion ni atipicos),
                suavizada con mediana movil e interpolada linealmente a traves del calentamiento.
ΔF/F0 corr.   = (F / (F0 * c(t)) - 1) * 100   (queda en 0 en los frames basales por construccion).
D (por ROI)   = mediana de ΔF/F0 corr. a T >= T_alta.
"""
import os
import numpy as np
import pandas as pd
from . import signal as sg


class Registro:
    """Un registro ya sincronizado: tiempo, temperatura y fluorescencia por frame, mascaras de calidad."""

    def __init__(self, fila, raiz, cfg):
        p = lambda c: os.path.join(raiz, fila[c])
        F = sg.leer_roi(p("roi_csv"))
        t, T, v = sg.leer_abf(p("abf"), cfg["ttl"]["canal_temperatura"], cfg["ttl"]["canal_ttl"])
        self.t, self.T = sg.temperatura_por_frame(t, T, v, fila["umbral_ttl"], len(F), cfg["ttl"].get("alineacion", "inicio"))
        self.n = len(self.T); self.F = F[: self.n]
        self.transicion = sg.leer_transicion(p("deriva_csv"), self.n)
        self.genotipo, self.registro = fila["genotipo"], int(fila["registro"])
        self.cfg = cfg
        nf = cfg["normalizacion"]["frames_f0"]
        self.F0 = self.F[:nf].mean(0); self.R = self.F / self.F0
        q = cfg["qc"]
        med = np.median(self.R, axis=1)
        res = med - pd.Series(med).rolling(q["ventana_atipicos"], center=True, min_periods=1).median().values
        sigma = 1.4826 * np.median(np.abs(res - np.median(res)))
        self.atipico = np.abs(res) > q["sigma_atipicos"] * sigma
        self.keep = ~self.transicion & ~self.atipico
        self.mediana = med
        self.T0 = float(np.median(self.T[:nf]))


def corregir(reg):
    """Devuelve (Fc [frames x ROIs] en %, c(t), mascara basal) o None si hay < 20 frames basales."""
    nrm = reg.cfg["normalizacion"]
    basal = reg.keep & (reg.T <= reg.T0 + nrm["margen_basal_C"])
    if basal.sum() < 20:
        return None
    idx = np.arange(reg.n)
    suave = pd.Series(reg.mediana[basal]).rolling(nrm["suavizado_basal"], center=True, min_periods=1).median().values
    c = np.interp(idx, idx[basal], suave)
    return (reg.R / c[:, None] - 1) * 100, c, basal


def metricas(reg):
    """Tablas por registro y por ROI + matriz ROI x bin de temperatura. Devuelve None si el registro no sirve."""
    cor = corregir(reg)
    if cor is None:
        return None
    Fc, c, basal = cor
    cfg = reg.cfg; g = cfg["graficos"]; nrm = cfg["normalizacion"]
    T, n = reg.T, reg.n; idx = np.arange(n)
    caliente = reg.keep & (T >= nrm["T_alta"])
    D = np.median(Fc[caliente], axis=0) if caliente.sum() >= 3 else np.full(Fc.shape[1], np.nan)
    ib = np.where(basal)[0]; tras = ib[ib > np.argmax(T > reg.T0 + nrm["margen_basal_C"])][-20:]
    post = Fc[tras].mean(0) if len(tras) >= 5 else np.full(Fc.shape[1], np.nan)
    inicio = int(np.argmax(T > reg.T0 + cfg["subida_previa"]["margen_inicio_C"]))
    pre = np.where(reg.keep & (idx < inicio))[0][-120:]
    pend = np.polyfit(reg.t[pre], reg.R[pre] * 100, 1)[0] if len(pre) >= 20 else np.full(Fc.shape[1], np.nan)
    bins = np.arange(g["bins_T"][0], g["bins_T"][1] + 1e-9, g["bins_T"][2]); centros = 0.5 * (bins[:-1] + bins[1:])
    ibin = np.digitize(T, bins) - 1
    M = np.full((Fc.shape[1], len(centros)), np.nan); cnt = np.zeros(len(centros), int)
    for b in range(len(centros)):
        m = reg.keep & (ibin == b); cnt[b] = m.sum()
        if m.sum():
            M[:, b] = Fc[m].mean(0)
    info = dict(genotipo=reg.genotipo, registro=reg.registro, rois=Fc.shape[1], frames=n, T0=round(reg.T0, 1), Tmax=round(float(T.max()), 1),
                frames_calientes=int(caliente.sum()), frames_transicion=int(reg.transicion.sum()), frames_atipicos=int(reg.atipico.sum()),
                deriva_comun_pct=round(float((c[-1] - 1) * 100), 1), D_mediana=round(float(np.nanmedian(D)), 2) if caliente.sum() >= 3 else np.nan)
    roi = pd.DataFrame(dict(genotipo=reg.genotipo, registro=reg.registro, roi=np.arange(1, Fc.shape[1] + 1), D=D,
                            ruido=Fc[basal].std(0), F0=reg.F0, post=post, pendiente_previa=pend))
    return info, roi, M, cnt, bins

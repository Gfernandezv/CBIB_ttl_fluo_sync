"""Promedio disparado por cada subida de temperatura: ¿hay un transitorio ligado al calentamiento?
Inicio = frame de maxima pendiente de la T suavizada dentro de cada subida (>= dT en `ventana` s), menos ~1 s.
A (por ROI) = media de ΔF/F_previo en [a, b] s tras el inicio. Nulo: la misma A con el mismo numero de inicios al azar.
p = fraccion del nulo >= A observada (prueba de aumento). Con azar puro ~alfa de los ROIs sale significativo."""
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter


def detectar_subidas(T, dt, dT, ventana, refractario):
    k = max(1, int(round(ventana / dt)))
    Ts = median_filter(T, size=3, mode="nearest")
    subida = np.zeros(len(T)); subida[: len(T) - k] = Ts[k:] - Ts[: len(T) - k]
    g = np.gradient(Ts); cand = np.where(subida >= dT)[0]
    inicios, i = [], 0
    while i < len(cand):
        i0 = cand[i]; seg = np.arange(i0, min(i0 + k + 1, len(T) - 1))
        inicios.append(max(0, seg[np.argmax(g[seg])] - int(round(1.0 / dt))))
        i = int(np.searchsorted(cand, i0 + k + int(round(refractario / dt))))
    return np.array(inicios)


def _alinear(F, on, npre, npost):
    out = np.empty((len(on), npre + npost + 1, F.shape[1]))
    for r, o in enumerate(on):
        out[r] = F[o - npre: o + npost + 1] / np.nanmean(F[o - npre: o], axis=0) - 1
    return out


def promedio_disparado(reg, cfg):
    d = cfg["disparado"]
    F = reg.F.copy(); F[reg.transicion] = np.nan          # solo frames de transicion (como en la v1)
    dt = float(np.median(np.diff(reg.t))); npre, npost = int(round(d["pre_s"] / dt)), int(round(d["post_s"] / dt))
    tt = np.arange(-npre, npost + 1) * dt; a, b = d["medir_s"]; ven = (tt >= a) & (tt <= b)
    inicios = detectar_subidas(reg.T, dt, d["dT"], d["ventana_subida_s"], d["refractario_s"])
    inicios = inicios[(inicios >= npre) & (inicios + npost < reg.n)]
    if len(inicios) < 2:
        return dict(subidas=len(inicios), A=np.full(F.shape[1], np.nan), p=np.full(F.shape[1], np.nan))
    A = np.nanmean(np.nanmean(_alinear(F, inicios, npre, npost)[:, ven, :], axis=1), axis=0)
    rng = np.random.default_rng(d["semilla"]); pool = np.arange(npre, reg.n - npost - 1)
    nulo = np.empty((d["n_nulo"], F.shape[1]))
    for i in range(d["n_nulo"]):
        on = np.sort(rng.choice(pool, len(inicios), replace=False))
        nulo[i] = np.nanmean(np.nanmean(_alinear(F, on, npre, npost)[:, ven, :], axis=1), axis=0)
    p = np.maximum((nulo >= A).mean(0), 1 / d["n_nulo"])
    return dict(subidas=len(inicios), A=A, p=p)

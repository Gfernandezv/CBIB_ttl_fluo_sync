"""Forma de la subida de F antes de calentar (T constante): recta vs exponencial saturante y = 1 + A(1 - exp(-t/tau)).
Mismos 2 parametros -> se compara el error cuadratico. Saturante = carga/iluminacion; lineal y sostenida = deriva/actividad."""
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import spearmanr


def _expo(t, A, tau):
    return 1 + A * (1 - np.exp(-t / tau))


def _ajustar(t, y):
    p = np.polyfit(t, y, 1); rl = float(np.sum((y - np.polyval(p, t)) ** 2))
    try:
        q, _ = curve_fit(_expo, t, y, p0=[0.05, 60.0], bounds=([-0.5, 3.0], [0.5, 3000.0]), maxfev=4000)
        re = float(np.sum((y - _expo(t, *q)) ** 2))
    except Exception:
        q, re = (np.nan, np.nan), np.nan
    return p, rl, q, re


def forma_subida(reg, cfg):
    c = cfg["subida_previa"]
    inicio = int(np.argmax(reg.T > reg.T0 + c["margen_inicio_C"]))
    pre = np.where(~reg.transicion & (np.arange(reg.n) < inicio))[0]
    if len(pre) < c["min_frames"]:
        return None, None
    tt = reg.t[pre] - reg.t[pre][0]; m = np.median(reg.R[pre], axis=1)
    p, rl, q, re = _ajustar(tt, m)
    mejora, amp = [], []
    for j in range(reg.R.shape[1]):
        pj, rlj, _, rej = _ajustar(tt, reg.R[pre, j]); amp.append(pj[0] * 100); mejora.append(rlj / rej if rej and rej > 0 else np.nan)
    dur = float(tt[-1])
    fila = dict(genotipo=reg.genotipo, registro=reg.registro, rois=reg.R.shape[1], s_previos=round(dur), pendiente_pct_s=round(p[0] * 100, 3),
                cambio_total_pct=round(p[0] * 100 * dur, 1), A_exp_pct=round(q[0] * 100, 1) if q[0] == q[0] else np.nan,
                tau_s=round(q[1]) if q[1] == q[1] else np.nan, RSS_recta_sobre_exp=round(rl / re, 2) if re and re == re else np.nan,
                pct_ROIs_exp_mejor=round(100 * float(np.nanmean(np.array(mejora) > 1.10)), 0),
                rho_pendiente_F0=round(float(spearmanr(np.array(amp), reg.F0)[0]), 2))
    return fila, (tt, m, p, q)


def curva(tt, p, q):
    return np.polyval(p, tt), (_expo(tt, *q) if q[0] == q[0] else None)

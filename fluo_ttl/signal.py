"""Lectura de un registro: abf (temperatura + TTL) y tabla de ROIs de ImageJ -> frames sincronizados."""
import numpy as np
import pandas as pd
import pyabf


def leer_abf(ruta, canal_T=0, canal_ttl=1):
    a = pyabf.ABF(str(ruta))
    a.setSweep(0, channel=canal_T); t = a.sweepX.copy(); T = a.sweepY.copy()
    a.setSweep(0, channel=canal_ttl); v = a.sweepY.copy()
    return t, T, v


def pulsos(v, umbral):
    """Flancos de subida/bajada del TTL por encima de `umbral`. Devuelve (inicios, finales)."""
    d = np.diff((v > umbral).astype(int))
    on, off = np.where(d == 1)[0] + 1, np.where(d == -1)[0] + 1
    n = min(len(on), len(off))
    return on[:n], off[:n]


def buscar_umbral(v, n_frames, n_cand=60):
    """Umbral TTL que da exactamente n_frames (o n_frames + 1) pulsos; el del medio de la meseta valida.
    Devuelve None si ninguno es consistente."""
    cand = np.linspace(np.percentile(v, 50) + 0.02, v.max() * 0.95, n_cand)
    ok = [c for c in cand if len(pulsos(v, c)[0]) in (n_frames, n_frames + 1)]
    return float(np.median(ok)) if ok else None


def temperatura_por_frame(t, T, v, umbral, n_frames):
    """Frame k = pulso TTL k; T del frame = promedio de la temperatura durante el pulso."""
    on, off = pulsos(v, umbral)
    m = min(len(on), n_frames)
    temp = np.array([T[o:f + 1].mean() for o, f in zip(on[:m], off[:m])])
    return t[on[:m]], temp


def leer_roi(ruta):
    df = pd.read_csv(ruta)
    cols = [c for c in df.columns if c.startswith("IntDen")]
    return df[cols].to_numpy(float)


def leer_transicion(ruta, n):
    return pd.read_csv(ruta)["transition"].to_numpy().astype(bool)[:n]

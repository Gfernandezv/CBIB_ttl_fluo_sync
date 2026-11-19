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


def temperatura_por_frame(t, T, v, umbral, n_frames, alineacion="inicio"):
    """Frame k = pulso TTL k; T del frame = promedio de la temperatura durante el pulso.
    Si hay mas pulsos que frames, `alineacion` elige cuales: 'inicio' = los primeros n, 'final' = los ultimos n."""
    on, off = pulsos(v, umbral)
    m = min(len(on), n_frames)
    if alineacion == "final":
        on, off = on[len(on) - m:], off[len(off) - m:]
    else:
        on, off = on[:m], off[:m]
    temp = np.array([T[o:f + 1].mean() for o, f in zip(on, off)])
    return t[on], temp


def leer_roi(ruta):
    df = pd.read_csv(ruta)
    cols = [c for c in df.columns if c.startswith("IntDen")]
    return df[cols].to_numpy(float)


def leer_transicion(ruta, n):
    return pd.read_csv(ruta)["transition"].to_numpy().astype(bool)[:n]


def diagnostico_pulsos(t, v, umbral, n_frames, T=None):
    """Describe los pulsos TTL frente a los frames para decidir por que no coinciden.
    Compara cada intervalo entre pulsos con el mediano y cada duracion con la mediana."""
    on, off = pulsos(v, umbral)
    n = len(on)
    out = dict(pulsos=n, frames=int(n_frames), dif=n - int(n_frames), razon=round(n / max(n_frames, 1), 3))
    if n < 3:
        return {**out, "sugerencia": "casi sin pulsos: revisa el umbral o el canal TTL"}
    dt = np.diff(t[on]); med = float(np.median(dt)); w = (t[off] - t[on])
    rel = dt / med; raros = np.where((rel < 0.7) | (rel > 1.3))[0]
    out.update(intervalo_mediano_s=round(med, 4), duracion_mediana_s=round(float(np.median(w)), 4),
               primer_intervalo=round(float(rel[0]), 2), ultimo_intervalo=round(float(rel[-1]), 2),
               n_intervalos_raros=int(len(raros)), posiciones_raras=",".join(str(i) for i in raros[:8]),
               pulsos_cortos=int((w < 0.5 * np.median(w)).sum()))
    d = out["dif"]
    if abs(out["razon"] - 2) < 0.1:
        s = "pulsos ~ 2x frames: cada frame puede generar 2 pulsos, o el stack tiene la mitad de los frames"
    elif d > 1 and rel[0] < 0.7 and len(raros) <= 2:
        s = "pulso(s) extra al INICIO (intervalo inicial corto)"
    elif d > 1 and rel[-1] < 0.7 and len(raros) <= 2:
        s = "pulso(s) extra al FINAL (intervalo final corto)"
    elif d > 1 and len(raros) == 0:
        s = "cadencia regular: el abf tiene mas pulsos que el stack (stack recortado o grabacion mas larga)"
    elif d > 1 and len(raros) > 2:
        s = "intervalos irregulares en medio: posibles pulsos espurios o pausas"
    elif d < -1:
        s = "faltan pulsos: el umbral puede ser demasiado alto"
    else:
        s = "dentro de la tolerancia"
    out["sugerencia"] = s
    if T is not None and 0 < d <= 20:
        Ta = np.array([T[o:f + 1].mean() for o, f in zip(on, off)])
        dd = np.abs(Ta[:n_frames] - Ta[d:])        # cambio de T de cada frame segun se alinee con los primeros o los ultimos pulsos
        out["sens_mediana_C"] = round(float(np.median(dd)), 2); out["sens_p95_C"] = round(float(np.percentile(dd, 95)), 2)
    return out

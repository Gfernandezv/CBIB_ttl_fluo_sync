"""Imagenes de microscopia: mapa espacial de ROIs coloreados por su respuesta (D) y vistas por temperatura.

Necesita, por registro, lo que dejan los macros de Fiji (todo es opcional; lo que falte se avisa):
  limpieza_ij/<nombre>_clean_view.tif | _clean_quant.tif   imagen (stack) alineada y recortada
  rois_ij/<nombre>_RoiSet.zip | RoiSet.zip                 contornos (mismo orden que las columnas de ROI_Results)
  rois_ij/<nombre>_rois_table_ij.csv | rois_table_ij.csv   centroide y area (alternativa si no hay RoiSet)
Los ROI estan en las coordenadas del stack alineado (clean_quant), no en las del original.
"""
import glob
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, Circle

from . import normalize


def _uno(carpeta, *patrones):
    for p in patrones:
        h = sorted(glob.glob(os.path.join(carpeta, p)))
        if h:
            return h[0]
    return None


def archivos_registro(fila, raiz):
    """Rutas de imagen/ROI de un registro (fila de registros.csv). None si no existen."""
    c = os.path.join(raiz, fila["carpeta"]); L, R = os.path.join(c, "limpieza_ij"), os.path.join(c, "rois_ij")
    return dict(stack=_uno(L, "*_clean_quant.tif", "*_clean_quant.tiff"), vista=_uno(L, "*_clean_view.tif", "*_clean_view.tiff"),
                roiset=_uno(R, "*_RoiSet.zip", "RoiSet.zip"), tabla=_uno(R, "*_rois_table_ij.csv", "rois_table_ij.csv"))


def leer_stack(ruta):
    import tifffile
    a = tifffile.imread(ruta)
    return a[None] if a.ndim == 2 else a


def leer_rois(arch):
    """Lista de dicts {x, y, poligono|None} en el orden de los ROI. Usa RoiSet.zip y, si no hay, la tabla de centroides."""
    if arch["roiset"]:
        import roifile
        rois = roifile.roiread(arch["roiset"]); rois = [rois] if isinstance(rois, roifile.ImagejRoi) else list(rois)
        out = []
        for r in rois:
            try:
                co = r.coordinates()
            except Exception:
                co = None
            if co is not None and len(co) >= 3:
                out.append(dict(x=float(co[:, 0].mean()), y=float(co[:, 1].mean()), poligono=np.asarray(co, float)))
            else:
                out.append(dict(x=(r.left + r.right) / 2, y=(r.top + r.bottom) / 2, poligono=None))
        return out
    if arch["tabla"]:
        t = pd.read_csv(arch["tabla"])
        return [dict(x=float(x), y=float(y), poligono=None) for x, y in zip(t["X"], t["Y"])]
    return None


def _imagen_fondo(arch, n=None):
    if arch["vista"]:
        v = leer_stack(arch["vista"]); return v.mean(0) if v.ndim == 3 and v.shape[0] > 1 else v[0]
    if arch["stack"]:
        s = leer_stack(arch["stack"]); return s[: (n or len(s))].mean(0)
    return None


def _norm(img, lo=1, hi=99.5):
    a, b = np.percentile(img, [lo, hi]); return np.clip((img - a) / max(b - a, 1e-9), 0, 1)


def _dibujar_rois(ax, rois, valores, vmax, resaltar=None, radio=5):
    cmap = plt.get_cmap("PuOr_r"); norm = matplotlib.colors.Normalize(-vmax, vmax); parches, col = [], []
    for r, v in zip(rois, valores):
        p = Polygon(r["poligono"], closed=True) if r["poligono"] is not None else Circle((r["x"], r["y"]), radio)
        parches.append(p); col.append(np.nan if v != v else v)
    pc = PatchCollection(parches, cmap=cmap, norm=norm, alpha=.55, edgecolor="k", linewidth=.3); pc.set_array(np.array(col)); ax.add_collection(pc)
    if resaltar is not None:
        pr = PatchCollection([p for p, h in zip(parches, resaltar) if h], facecolor="none", edgecolor="red", linewidth=1.3); ax.add_collection(pr)
    return pc


def mapa_espacial(fila, raiz, por_roi, ax=None, vmax=5.0, p99=None, titulo=None):
    """Imagen del registro con cada ROI coloreado por D (naranja = aumenta, violeta = disminuye).
    Si se da `p99`, los ROI con D > p99 llevan contorno rojo. Devuelve (ax, motivo) y motivo = None si se pudo dibujar."""
    arch = archivos_registro(fila, raiz); rois = leer_rois(arch)
    if rois is None:
        return ax, "faltan RoiSet.zip y rois_table_ij.csv en rois_ij/"
    d = por_roi[(por_roi.genotipo == fila["genotipo"]) & (por_roi.registro == fila["registro"])].sort_values("roi").D.values
    if len(d) != len(rois):
        return ax, f"hay {len(rois)} ROI en el RoiSet/tabla y {len(d)} en por_roi.csv: no se puede emparejar"
    img = _imagen_fondo(arch)
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))
    if img is not None:
        ax.imshow(_norm(img), cmap="gray", interpolation="nearest")
    else:
        xs = [r["x"] for r in rois]; ys = [r["y"] for r in rois]; ax.set_xlim(min(xs) - 20, max(xs) + 20); ax.set_ylim(max(ys) + 20, min(ys) - 20); ax.set_facecolor("0.1")
    pc = _dibujar_rois(ax, rois, d, vmax, None if p99 is None else d > p99)
    ax.set_title(titulo or f"{fila['genotipo']} reg {fila['registro']} — D por ROI", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    plt.colorbar(pc, ax=ax, fraction=.046, pad=.02, label="D = ΔF/F₀ corr. a T alta (%)")
    return ax, None


def top_rois(por_roi, genotipo, registro, n=6, modo="sube"):
    """Numeros de ROI (1..N) con mayor (modo='sube') o menor ('baja') D en un registro."""
    d = por_roi[(por_roi.genotipo == genotipo) & (por_roi.registro == registro)].dropna(subset=["D"])
    return d.sort_values("D", ascending=(modo == "baja")).roi.head(n).tolist()


def vistas_por_temperatura(fila, raiz, cfg, rois_n, temps=(22, 26, 30, 34, 38), ancho_bin=1.0, caja=40, ruta=None):
    """Para cada ROI (numero 1..N) y cada temperatura: recorte de la media de los frames con T dentro de ±ancho_bin/2.
    Cada fila usa su propia escala de gris (comparable entre columnas, no entre filas)."""
    arch = archivos_registro(fila, raiz); rois = leer_rois(arch)
    if not arch["stack"] or rois is None:
        return None, "faltan el stack *_clean_quant.tif o los ROI"
    reg = normalize.Registro(fila, raiz, cfg); st = leer_stack(arch["stack"])
    if st.shape[0] < reg.n:
        return None, f"el stack tiene {st.shape[0]} frames y el registro {reg.n}"
    st = st[: reg.n]; h = caja // 2
    fig, axs = plt.subplots(len(rois_n), len(temps), figsize=(1.9 * len(temps), 1.9 * len(rois_n)), squeeze=False)
    for i, n in enumerate(rois_n):
        r = rois[n - 1]; cx, cy = int(round(r["x"])), int(round(r["y"])); y0, x0 = max(cy - h, 0), max(cx - h, 0)
        crops = []
        for T in temps:
            m = reg.keep & (np.abs(reg.T - T) <= ancho_bin / 2)
            crops.append(st[m][:, y0:cy + h, x0:cx + h].mean(0) if m.sum() >= 3 else None)
        valid = [c for c in crops if c is not None]
        lo, hi = (np.percentile(np.concatenate([c.ravel() for c in valid]), [1, 99.5]) if valid else (0, 1))
        for j, (T, c) in enumerate(zip(temps, crops)):
            ax = axs[i, j]; ax.set_xticks([]); ax.set_yticks([])
            if c is None:
                ax.set_facecolor("0.85"); ax.text(.5, .5, "sin frames", ha="center", va="center", fontsize=7, transform=ax.transAxes)
            else:
                ax.imshow(c, cmap="gray", vmin=lo, vmax=hi, interpolation="nearest")
            if i == 0:
                ax.set_title(f"{T:g} °C", fontsize=9)
            if j == 0:
                ax.set_ylabel(f"ROI {n}", fontsize=9)
    fig.suptitle(f"{fila['genotipo']} reg {fila['registro']}: media de frames por temperatura", fontsize=10); fig.tight_layout()
    if ruta:
        fig.savefig(ruta, dpi=130, bbox_inches="tight")
    return fig, None

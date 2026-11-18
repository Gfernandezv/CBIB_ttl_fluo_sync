"""Figuras estandar: pseudocolor ROI x temperatura (agrupado y por muestra), comparacion entre genotipos, forma de la subida previa."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from . import preheat


def _cmap():
    c = plt.get_cmap("PuOr_r").copy(); c.set_bad("0.85"); return c


def _ordenar(M, cnt, centros, minf):
    Mm = M.copy(); Mm[:, cnt < minf] = np.nan
    alta = np.nanmean(Mm[:, centros >= 30], axis=1)
    return Mm, np.argsort(-np.nan_to_num(alta, nan=-99))        # mayor aumento arriba


def pseudocolor_agrupado(geno, regs, bins, cfg, ruta):
    g = cfg["graficos"]; vmax = g["escala_color_pct"]; centros = 0.5 * (bins[:-1] + bins[1:])
    M = np.vstack([r["M"] for r in regs.values()]); cnt = np.sum([r["cnt"] for r in regs.values()], axis=0)
    Mm, orden = _ordenar(M, cnt, centros, g["min_frames_bin"])
    fig = plt.figure(figsize=(9.5, 8.5)); gs = fig.add_gridspec(3, 2, height_ratios=[0.9, 3.6, 0.7], width_ratios=[40, 1.2], hspace=0.12, wspace=0.03)
    a0 = fig.add_subplot(gs[0, 0]); a1 = fig.add_subplot(gs[1, 0], sharex=a0); a2 = fig.add_subplot(gs[2, 0], sharex=a0)
    q1, med, q3 = np.nanpercentile(Mm, [25, 50, 75], axis=0)
    a0.fill_between(centros, q1, q3, color="0.8", step="mid"); a0.plot(centros, med, "k-o", ms=3); a0.axhline(0, c="gray", ls="--", lw=.8)
    a0.set_ylabel("mediana (IQR)\nΔF/F₀ corr. %"); a0.set_title(f"{geno}: {len(regs)} registros, {M.shape[0]} ROIs — ΔF/F₀ (%) corregido por deriva", fontsize=10)
    plt.setp(a0.get_xticklabels(), visible=False)
    im = a1.imshow(Mm[orden], aspect="auto", cmap=_cmap(), vmin=-vmax, vmax=vmax, interpolation="nearest", extent=[bins[0], bins[-1], Mm.shape[0], 0])
    a1.set_ylabel("ROIs ordenados por ΔF/F₀ medio a T ≥ 30 °C\n(mayor aumento arriba)"); plt.setp(a1.get_xticklabels(), visible=False)
    a2.bar(centros, cnt, width=1.0, color=np.where(cnt < g["min_frames_bin"], "C3", "0.4")); a2.set_yscale("log"); a2.set_ylim(0.8, None)
    a2.set_ylabel("frames\npor bin"); a2.set_xlabel("Temperatura (°C)")
    for a in (a0, a1, a2):
        a.axvline(g["linea_C"], c="k", ls="--", lw=1.2)
    a1.set_xlim(bins[0], bins[-1]); plt.colorbar(im, cax=fig.add_subplot(gs[1, 1]), label="ΔF/F₀ corr. (%)   [naranja = gana · violeta = pierde]")
    fig.savefig(ruta, dpi=130, bbox_inches="tight"); plt.close(fig)


def pseudocolor_por_muestra(geno, regs, bins, cfg, ruta):
    g = cfg["graficos"]; vmax = g["escala_color_pct"]; centros = 0.5 * (bins[:-1] + bins[1:])
    nm = len(regs); nc = min(nm, 4); nr = int(np.ceil(nm / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(4.4 * nc + 1.5, 4.2 * nr), sharex=True, squeeze=False)
    for ax in axs.ravel():
        ax.axis("off")
    for ax, (s, r) in zip(axs.ravel(), regs.items()):
        ax.axis("on"); Mm, o = _ordenar(r["M"], r["cnt"], centros, g["min_frames_bin"])
        im = ax.imshow(Mm[o], aspect="auto", cmap=_cmap(), vmin=-vmax, vmax=vmax, interpolation="nearest", extent=[bins[0], bins[-1], Mm.shape[0], 0])
        ax.axvline(g["linea_C"], c="k", ls="--", lw=1); ax.set_title(f"{geno} reg {s} (n={Mm.shape[0]})", fontsize=10); ax.set_xlim(bins[0], bins[-1])
    for ax in axs[-1]:
        ax.set_xlabel("Temperatura (°C)")
    for ax in axs[:, 0]:
        ax.set_ylabel("ROIs")
    fig.colorbar(im, ax=axs, label="ΔF/F₀ corr. (%)", shrink=.8, pad=.01)
    fig.suptitle(f"ΔF/F₀ corregido por deriva, {geno} por muestra (gris = < {g['min_frames_bin']} frames)")
    fig.savefig(ruta, dpi=110, bbox_inches="tight"); plt.close(fig)


def comparacion(res, rois, bins, cfg, ruta, p1, p99):
    g = cfg["graficos"]; centros = 0.5 * (bins[:-1] + bins[1:]); cols = plt.cm.tab10(np.arange(10)); ref = cfg["referencia"]
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.5))
    for k, (geno, regs) in enumerate(res.items()):
        M = np.vstack([r["M"] for r in regs.values()]); cnt = np.sum([r["cnt"] for r in regs.values()], axis=0); M[:, cnt < g["min_frames_bin"]] = np.nan
        ax[0].plot(centros, np.nanmedian(M, 0), "-o", ms=3, c=cols[k], lw=2.2 if geno == ref else 1.4, label=f"{geno} (n={M.shape[0]})")
        d = np.sort(rois[rois.genotipo == geno].D.dropna().values)
        ax[1].plot(d, np.arange(1, len(d) + 1) / max(len(d), 1), c=cols[k], lw=2.2 if geno == ref else 1.4, label=geno)
    ax[0].axhline(0, c="gray", ls="--", lw=.8); ax[0].axvline(g["linea_C"], c="k", ls="--", lw=1)
    ax[0].set_xlabel("Temperatura (°C)"); ax[0].set_ylabel("mediana ΔF/F₀ corregido (%)"); ax[0].legend(fontsize=8); ax[0].set_title("Mediana de ROIs vs temperatura, por genotipo")
    ax[1].axvline(p1, c="k", ls=":", lw=1); ax[1].axvline(p99, c="k", ls=":", lw=1, label=f"P1 / P99 de {ref}"); ax[1].set_xlim(-6, 8)
    ax[1].set_xlabel("D = ΔF/F₀ corregido mediano a T alta, por ROI (%)"); ax[1].set_ylabel("fracción acumulada de ROIs"); ax[1].legend(fontsize=8)
    ax[1].set_title("Distribución por ROI del efecto a T alta"); fig.tight_layout(); fig.savefig(ruta, dpi=120); plt.close(fig)


def forma_subida(sub, trazas, ruta):
    if len(sub) == 0:
        return
    nr = int(np.ceil(len(sub) / 4)); fig, axs = plt.subplots(nr, 4, figsize=(17, 3.4 * nr), squeeze=False)
    for ax in axs.ravel():
        ax.axis("off")
    for ax, (_, r) in zip(axs.ravel(), sub.iterrows()):
        tt, m, p, q = trazas[(r.genotipo, r.registro)]; ax.axis("on"); lin, exp = preheat.curva(tt, p, q)
        ax.plot(tt, (m - 1) * 100, "k", lw=.8, label="mediana ROIs"); ax.plot(tt, (lin - 1) * 100, "C0", label="recta")
        if exp is not None:
            ax.plot(tt, (exp - 1) * 100, "C3", label=f"exp τ={q[1]:.0f}s")
        ax.set_title(f"{r.genotipo} reg {r.registro} (RSS recta/exp = {r.RSS_recta_sobre_exp})", fontsize=9)
        ax.set_xlabel("s desde el inicio (T constante)"); ax.set_ylabel("F/F₀ − 1 (%)"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(ruta, dpi=110); plt.close(fig)

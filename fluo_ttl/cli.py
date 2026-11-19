"""Linea de comandos:  fluo-ttl scan | run | reproduce"""
import argparse
import warnings
import json
import os
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd
import yaml

from . import config as cfgmod, registry, normalize, triggered, preheat, plots, provenance

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLAS = ["registros.csv", "por_registro.csv", "por_roi.csv", "tabla_genotipos.csv", "disparado.csv", "forma_subida.csv"]


def _scan(raiz, cfg, out):
    os.makedirs(out, exist_ok=True)
    reg = registry.escanear(raiz, cfg)
    reg.to_csv(os.path.join(out, "registros.csv"), index=False)
    return reg


def _sin_registros(reg, raiz, cfg):
    """Mensaje claro cuando el escaneo no deja nada que analizar (en vez de fallar mas adelante)."""
    if len(reg) == 0:
        sub = sorted(d for d in os.listdir(raiz) if os.path.isdir(os.path.join(raiz, d))) if os.path.isdir(raiz) else []
        extra = [d for d in sub if d not in cfg["genotipos"]]
        msg = [f"No se encontro ningun registro en: {raiz}",
               "Se busca  <datos>/<genotipo>/**/rois_ij/*_ROI_Results.csv  junto a un .abf y limpieza_ij/*drift_shifts_ij.csv.",
               f"Genotipos buscados (configs): {', '.join(cfg['genotipos'])}"]
        if not os.path.isdir(raiz):
            msg.append("La carpeta no existe: revisa la ruta (entre comillas si tiene espacios).")
        elif extra:
            msg.append("Carpetas dentro de la ruta que NO son genotipos conocidos: " + ", ".join(extra[:15]))
            msg.append("Si una de ellas contiene los genotipos, usa esa como carpeta de datos; si son genotipos nuevos, agregalos con 'genotipos:' en un YAML (-c).")
        else:
            msg.append("La carpeta no contiene subcarpetas de genotipo.")
        sys.exit("\n".join(msg))
    if not reg.incluido.any():
        sys.exit("Hay %d registros pero ninguno entra al analisis:\n%s" % (len(reg), reg.groupby("motivo").size().to_string()))


def cmd_scan(a):
    cfg = cfgmod.cargar(a.config); reg = _scan(a.datos, cfg, a.salida)
    if len(reg) == 0:
        _sin_registros(reg, a.datos, cfg)
    print(reg[["genotipo", "registro", "frames", "pulsos", "umbral_ttl", "incluido", "motivo"]].to_string(index=False))
    print(f"\n{reg.incluido.sum()} incluidos, {(~reg.incluido).sum()} excluidos -> {a.salida}/registros.csv")
    if a.detalle:
        det = registry.detalle_pulsos(reg, a.datos, cfg); det.to_csv(os.path.join(a.salida, "pulsos_detalle.csv"), index=False)
        pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
        print("\nDiagnostico de pulsos TTL (registros con diferencia > 1):")
        x = det[det.dif.abs() > 1]
        print(x[["genotipo", "registro", "pulsos", "frames", "dif", "primer_intervalo", "ultimo_intervalo", "n_intervalos_raros", "sugerencia"]].to_string(index=False) if len(x) else "  todos coinciden")
        print(f"-> {a.salida}/pulsos_detalle.csv")


def ejecutar(raiz, cfg, out, repo_dir=REPO):
    """Corre todo el analisis; deja tablas, figuras y provenance.json en `out`."""
    reg = _scan(raiz, cfg, out)
    _sin_registros(reg, raiz, cfg)
    figs = os.path.join(out, "figuras"); os.makedirs(figs, exist_ok=True)
    infos, rois, trig, subs, res, trazas = [], [], [], [], {}, {}
    bins = None
    for _, f in reg[reg.incluido].iterrows():
        r = normalize.Registro(f, raiz, cfg); m = normalize.metricas(r)
        if m is None:
            print(f"  {f.genotipo} {f.registro}: < 20 frames basales, excluido"); continue
        info, roi, M, cnt, bins = m
        t = triggered.promedio_disparado(r, cfg)
        info.update(subidas_calor=t["subidas"]); infos.append(info)
        roi = roi.assign(A_disparado=t["A"], p_disparado=t["p"]); rois.append(roi)
        res.setdefault(f.genotipo, {})[f.registro] = dict(M=M, cnt=cnt)
        fila, tr = preheat.forma_subida(r, cfg)
        if fila:
            subs.append(fila); trazas[(f.genotipo, f.registro)] = tr
        print(f"  {f.genotipo} {f.registro}: {info['rois']} ROIs, D mediana {info['D_mediana']}")
    por_reg = pd.DataFrame(infos); por_roi = pd.concat(rois, ignore_index=True); sub = pd.DataFrame(subs)
    ref = cfg["referencia"]; d_ref = por_roi[por_roi.genotipo == ref].D.dropna()
    if len(d_ref) < 20:
        raise SystemExit(f"Hay muy pocos ROIs de la referencia '{ref}' para fijar umbrales.")
    p1, p99 = np.percentile(d_ref, [1, 99])
    filas = []
    for g, d in por_roi.groupby("genotipo"):
        x = d.D.dropna()
        filas.append(dict(genotipo=g, registros=d.registro.nunique(), rois=len(x), D_mediana=round(x.median(), 2),
                          pct_sobre_P99=round(100 * (x > p99).mean(), 1), pct_bajo_P1=round(100 * (x < p1).mean(), 1),
                          P1_ref=round(p1, 2), P99_ref=round(p99, 2)))
    tab = pd.DataFrame(filas)
    trig = por_roi[["genotipo", "registro", "roi", "A_disparado", "p_disparado"]]
    por_reg.to_csv(os.path.join(out, "por_registro.csv"), index=False); por_roi.drop(columns=["A_disparado", "p_disparado"]).to_csv(os.path.join(out, "por_roi.csv"), index=False)
    tab.to_csv(os.path.join(out, "tabla_genotipos.csv"), index=False); trig.to_csv(os.path.join(out, "disparado.csv"), index=False)
    sub.to_csv(os.path.join(out, "forma_subida.csv"), index=False)
    for g, regs in res.items():
        plots.pseudocolor_agrupado(g, regs, bins, cfg, os.path.join(figs, f"pseudocolor_{g}_agrupado.png"))
        plots.pseudocolor_por_muestra(g, regs, bins, cfg, os.path.join(figs, f"pseudocolor_{g}_porMuestra.png"))
    plots.comparacion(res, por_roi, bins, cfg, os.path.join(figs, "comparacion_genotipos.png"), p1, p99)
    plots.forma_subida(sub, trazas, os.path.join(figs, "forma_subida.png"))
    provenance.guardar(provenance.construir(raiz, reg, cfg, repo_dir), os.path.join(out, "provenance.json"))
    with open(os.path.join(out, "config_usada.yaml"), "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)
    print("\n" + tab.to_string(index=False)); print(f"\nListo -> {out}")
    return tab


def cmd_run(a):
    ejecutar(a.datos, cfgmod.cargar(a.config), a.salida)


def cmd_reproduce(a):
    """Verifica hashes de entrada, vuelve a correr con la config guardada y compara las tablas con las originales."""
    prov = json.load(open(os.path.join(a.salida, "provenance.json"), encoding="utf-8"))
    mal = [r for r, h in prov["entradas_sha256"].items() if not os.path.exists(os.path.join(a.datos, r)) or provenance.sha256(os.path.join(a.datos, r)) != h]
    if mal:
        print("Los datos de entrada no coinciden con los de la corrida original:"); [print("  ", m) for m in mal[:20]]; sys.exit(2)
    print(f"Entradas identicas ({len(prov['entradas_sha256'])} archivos). Codigo original: {prov['commit']} (fluo_ttl {prov['fluo_ttl']})")
    tmp = tempfile.mkdtemp(prefix="fluo_ttl_repro_")
    try:
        ejecutar(a.datos, prov["config"], tmp)
        ok = True
        for t in TABLAS:
            x, y = pd.read_csv(os.path.join(a.salida, t)), pd.read_csv(os.path.join(tmp, t))
            try:
                pd.testing.assert_frame_equal(x, y, check_exact=False, rtol=1e-9, atol=1e-9); print(f"  {t}: identica")
            except AssertionError as e:
                ok = False; print(f"  {t}: DIFIERE\n     {str(e)[:300]}")
        print("\nREPRODUCIBLE" if ok else "\nNO reproduce"); sys.exit(0 if ok else 1)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main(argv=None):
    warnings.filterwarnings("ignore", message="All-NaN slice", category=RuntimeWarning)
    p = argparse.ArgumentParser(prog="fluo-ttl", description="Fluorescencia (Fluo-4) sincronizada con temperatura (abf/TTL)")
    s = p.add_subparsers(dest="cmd", required=True)
    for nombre, fn, ayuda in (("scan", cmd_scan, "lista registros, umbral TTL y exclusiones"), ("run", cmd_run, "analisis completo"), ("reproduce", cmd_reproduce, "reejecuta y compara con una corrida previa")):
        q = s.add_parser(nombre, help=ayuda); q.set_defaults(fn=fn)
        q.add_argument("datos", help="carpeta raiz: <datos>/<genotipo>/**/{abf, rois_ij, limpieza_ij}")
        q.add_argument("-o", "--salida", default="resultados", help="carpeta de salida (default: resultados)")
        if nombre != "reproduce":
            q.add_argument("-c", "--config", help="YAML con parametros que sobreescriben fluo_ttl/default.yaml")
        if nombre == "scan":
            q.add_argument("--detalle", action="store_true", help="diagnostico de pulsos TTL (intervalos, pulsos extra al inicio/final)")
    a = p.parse_args(argv); a.fn(a)


if __name__ == "__main__":
    main()

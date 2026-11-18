"""Procedencia de una corrida: version del codigo (commit), versiones de paquetes, config y hash de cada archivo de entrada."""
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata

from . import __version__

PAQUETES = ["numpy", "pandas", "scipy", "matplotlib", "pyabf", "pyyaml"]


def sha256(ruta, bloque=1 << 20):
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for b in iter(lambda: f.read(bloque), b""):
            h.update(b)
    return h.hexdigest()


def _git(cmd, cwd):
    try:
        return subprocess.run(["git"] + cmd, cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


def construir(raiz, registro, cfg, repo_dir):
    entradas = {}
    for _, f in registro[registro.incluido].iterrows():
        for col in ("abf", "roi_csv", "deriva_csv"):
            entradas[f[col]] = sha256(os.path.join(raiz, f[col]))
    return dict(
        fluo_ttl=__version__, commit=_git(["rev-parse", "HEAD"], repo_dir), cambios_sin_commit=bool(_git(["status", "--porcelain", "--", "fluo_ttl", "configs"], repo_dir)),
        fecha_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"), python=sys.version.split()[0], sistema=platform.platform(),
        paquetes={p: metadata.version(p) for p in PAQUETES}, config=cfg, entradas_sha256=entradas)


def guardar(prov, ruta):
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=2, ensure_ascii=False, sort_keys=True)

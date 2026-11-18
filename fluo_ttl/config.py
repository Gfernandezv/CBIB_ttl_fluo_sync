import copy
from pathlib import Path
import yaml

DEFECTO = Path(__file__).resolve().parent / "default.yaml"


def _fusionar(base, extra):
    for k, v in extra.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _fusionar(base[k], v)
        else:
            base[k] = v
    return base


def cargar(ruta=None):
    """Lee default.yaml y le superpone el YAML del usuario (solo las claves que cambia)."""
    cfg = yaml.safe_load(DEFECTO.read_text(encoding="utf-8"))
    if ruta:
        cfg = _fusionar(copy.deepcopy(cfg), yaml.safe_load(Path(ruta).read_text(encoding="utf-8")) or {})
    return cfg

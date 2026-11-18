# fluo_ttl 2.0

Respuesta de fluorescencia (Fluo-4) de TRPM3 sincronizada con temperatura. Cada fotograma de imagen se asigna a su pulso TTL en el registro `.abf`, de donde se toma la temperatura. Se calcula un ΔF/F₀ corregido por la deriva común de cada registro y se compara el efecto a temperatura alta entre genotipos. Una corrida deja tablas, figuras y un `provenance.json` que permite repetirla y comprobar que da lo mismo.

> La v1 (notebooks 01–04 y scripts) está en `legacy/`. Sigue funcionando pero normaliza distinto (ver «Qué cambió»).

## Instalación

```bash
conda env create -f environment.yml && conda activate ttl_imagej   # o: pip install -e .
fluo-ttl --help
```

## Organización de los datos

La carpeta de primer nivel es el **genotipo** (`Wt` es la referencia). Dentro, cada registro tiene su carpeta con tres cosas:

```text
<datos>/
├── Wt/Nov06/Nov062026_2/
│   ├── 26n06000.abf                              # canal 0 = temperatura, canal 1 = TTL
│   ├── rois_ij/Nov062026_2_Stack_ROI_Results.csv   # medición de macros/03_medicion.ijm
│   └── limpieza_ij/Nov062026_2_Stack_drift_shifts_ij.csv   # marca los frames de transición
├── 1236/…   1250/…   1257/…   1227/…   1243/…   1265/…
```

Los genotipos reconocidos y los demás parámetros están en `fluo_ttl/default.yaml`.

## Uso

```bash
fluo-ttl scan  <datos> -o resultados            # lista registros, umbral TTL, pulsos vs frames y exclusiones
fluo-ttl run   <datos> -o resultados            # análisis completo
fluo-ttl run   <datos> -o resultados -c mi_experimento.yaml   # parámetros propios (solo hace falta lo que cambia)
fluo-ttl reproduce <datos> -o resultados        # verifica hashes de entrada, reejecuta y compara tablas
```

`mi_experimento.yaml` (ver `configs/ejemplo_experimento.yaml`) se fusiona sobre `fluo_ttl/default.yaml`; por ejemplo:

```yaml
referencia: Wt
ttl: {umbral_fijo: {Wt: 0.1}}
normalizacion: {frames_f0: 30}
```

### Salidas (`resultados/`)

| Archivo | Contenido |
|---|---|
| `registros.csv` | un renglón por registro: frames, pulsos, umbral TTL, incluido y motivo de exclusión |
| `por_registro.csv`, `por_roi.csv` | deriva común, D por ROI, ruido basal, pendiente previa al calentamiento |
| `tabla_genotipos.csv` | ROIs, mediana de D y % por encima de P99 / por debajo de P1 de la referencia |
| `disparado.csv` | promedio disparado por los ascensos de temperatura (A y p frente a un nulo aleatorio) |
| `forma_subida.csv` | subida previa al calor: recta vs exponencial saturante |
| `figuras/` | pseudocolor agrupado y por muestra, comparación entre genotipos, forma de la subida |
| `provenance.json`, `config_usada.yaml` | commit, versiones, parámetros y SHA-256 de cada archivo de entrada |

## Método

1. **Sincronía.** Fotograma *k* = pulso TTL *k*; su temperatura es la media del canal 0 durante el pulso. El umbral TTL se busca automáticamente (o se fija por genotipo). Un registro se excluye si |pulsos − frames| > 1; la exclusión queda registrada con su motivo.
2. **Calidad.** Se descartan los frames de transición (tabla de deriva de ImageJ) y los atípicos globales (mediana de ROIs frente a mediana móvil de 11 frames, > 6σ).
3. **ΔF/F₀ corregido.** F₀ = media de los primeros 20 frames de cada ROI. La deriva común c(t) es la mediana entre ROIs de F/F₀ en frames basales (T ≤ T₀ + 2 °C), suavizada e interpolada a través del calentamiento. ΔF/F₀ = (F / (F₀·c(t)) − 1)·100; positivo (naranja) = aumento.
4. **Efecto a temperatura alta.** D = mediana de ΔF/F₀ corregido a T ≥ 36 °C por ROI. Los umbrales de «respuesta» son P1 y P99 de D en la referencia (Wt).
5. **Controles.** Promedio disparado por ascensos de temperatura con 1000 inicios aleatorios como nulo (semilla fija); forma de la subida previa al calor.

## Qué cambió respecto a la v1

La v1 calculaba F₀ como la intensidad promediada a la temperatura objetivo ±1 °C sobre **toda** la grabación, de modo que mezclaba frames tempranos y tardíos y la deriva (blanqueo o aumento lento) aparecía como respuesta. La v2 fija F₀ en el inicio y resta la deriva común. Las clases de la v1 no son comparables con las de la v2.

## Limitaciones conocidas

- Siempre hay un pulso TTL extra respecto de los frames (se tolera 1); no está resuelto si es el primero o el último.
- Los registros cuyos pulsos y frames no coinciden se excluyen (p. ej. Wt 13 y varios de 1265); el emparejamiento imagen–abf de 1265 debe revisarse.
- Ningún genotipo muestra aumento de fluorescencia sincronizado con el calor; varias grabaciones de mutantes suben antes de calentar, por lo que D de esos ROIs sigue su pendiente previa. Hace falta un periodo previo más largo a temperatura constante y un control positivo.

## Repositorio

`fluo_ttl/` paquete · `configs/` ejemplo de config de experimento · `macros/` Fiji (limpieza → ROIs → medición) · `notebooks/05_resultados_v2.ipynb` resultados · `legacy/` v1.

# download-genome

Herramienta en línea de comandos que, a partir de un TSV con genomas y accesiones GenBank, descarga de **NCBI** archivos **GenBank** (`.gb`) y **GFF3** (`.gff3`), con pausa entre peticiones y reporte opcional de errores.

![Python](https://img.shields.io/badge/python-3.11+-blue.svg)
![uv](https://img.shields.io/badge/uv-enabled-5c6cbc.svg)
![NCBI](https://img.shields.io/badge/data-NCBI_Entrez-0071bc.svg)

## Requisitos

- **Python** >= 3.11 (recomendado: la versión en `.python-version`).
- [uv](https://docs.astral.sh/uv/) para el entorno e instalación.
- Conexión a internet (NCBI Entrez y visor sviewer).
- Un correo electrónico válido: NCBI exige identificarte al usar Entrez; opcionalmente una **API key** de NCBI para mayor cuota.

### Dependencias de Python

| Paquete   | Uso en el proyecto         |
|----------|----------------------------|
| pandas   | Leer TSV, escribir reporte |
| requests | Descargar GFF3 vía HTTP   |
| biopython| `Bio.Entrez` (efetch, esearch) |

Están declaradas en `pyproject.toml` e instalables con `uv sync`.

## Instalación con uv

Clona el repositorio y entra al directorio del proyecto.

```bash
cd download_genome
```

Instala el proyecto y dependencias en un entorno virtual (crea/actualiza el venv según el lockfile):

```bash
uv sync
```

Si añades dependencias a mano en el futuro:

```bash
uv add "paquete>=x.y"
```

Ejecución con el intérprete del entorno de uv (sin activar el venv a mano):

```bash
uv run python main.py --help
```

## Formato de entrada (TSV)

El archivo de entrada es **separado por tabulaciones**, con al menos estas columnas (nombres por defecto; puedes cambiarlos con flags):

| Columna (por defecto) | Contenido |
|------------------------|-----------|
| `genome_id`            | Identificador del genoma |
| `genome_name`          | Nombre legible |
| `genbank_accessions`   | Una o más accesiones separadas por **comas** en la misma celda |

Si los nombres de columna difieren, usa `--col-genome-id`, `--col-genome-name` y `--col-accessions`.

## Uso

### Ayuda

```bash
uv run python main.py --help
```

### Ejemplo mínimo

NCBI requiere `--email` (usa un correo de contacto real).

```bash
uv run python main.py --email tu.correo@institucion.org -i genomas.tsv -o descargas_ncbi
```

### Parámetros útiles

- `-i` / `--input` — TSV (por defecto: `genomas.tsv`).
- `-o` / `--output-dir` — Directorio base de salida (por defecto: `descargas_ncbi`).
- `--delay` — Segundos entre peticiones (sin API key de NCBI, conviene ≥ **0.34**).
- `--api-key` — API key opcional de NCBI.
- Columnas: `--col-genome-id`, `--col-genome-name`, `--col-accessions`.

## Salida

Bajo el directorio de `--output-dir` se crean, entre otros:

- `genbank/` — archivos `{accession}.gb`
- `gff3/` — archivos `{accession}.gff3`
- `reporte_errores.tsv` — solo si hubo incidencias (accesion inexistente, fallos de descarga, etc.)

## Política de NCBI

Cumple las [directrices de NCBI para el uso de E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/): no abuses la frecuencia de peticiones; ajusta `delay` y, con uso intenso, registro e **API key**.

## Topics sugeridos (GitHub)

Añade en el repositorio: **Settings → General → Topics**.

`python` · `ncbi` · `entrez` · `genbank` · `gff3` · `bioinformatics` · `genomics` · `pandas` · `uv` · `regulondb`

## Licencia

Especifica la licencia de este repositorio (p. ej. en `LICENSE`) o actualiza este apartado.

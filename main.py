"""Descarga GenBank y GFF3 desde NCBI para genomas listados en un TSV.

Lee un archivo tabular con identificadores de genoma, nombres y accesiones
GenBank (separadas por comas). Para cada accession verifica existencia en
NCBI, descarga ``.gb`` y ``.gff3``, y opcionalmente escribe un reporte de
errores.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from Bio import Entrez


def resolve_column_label(spec: str, columns: pd.Index) -> str:
    """Obtiene la etiqueta de columna por nombre o por índice en base 1.

    Si ``spec`` coincide con un nombre de columna existente, se usa tal cual
    (así una columna llamada ``\"4\"`` no se confunde con el índice 4).
    Si es solo dígitos y no hay columna con ese nombre, se interpreta como
    número de columna 1-based (1 = primera columna).

    Args:
        spec: Nombre de cabecera o número (ej. ``\"4\"`` para la cuarta columna).
        columns: ``df.columns`` del DataFrame leído.

    Returns:
        Etiqueta real en ``columns`` para indexar filas.

    Raises:
        ValueError: Índice fuera de rango o columna inexistente.
    """
    name = spec.strip()
    if name in columns:
        return name
    if name.isdigit():
        idx = int(name)
        if 1 <= idx <= len(columns):
            return columns[idx - 1]
        raise ValueError(
            f"Índice de columna inválido {idx!r}: el archivo tiene "
            f"{len(columns)} columnas (usa 1..{len(columns)})."
        )
    raise ValueError(f"No existe la columna {name!r}. Cabeceras: {list(columns)}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Construye el parser CLI y devuelve los argumentos parseados."""
    parser = argparse.ArgumentParser(
        description=(
            "Descarga archivos GenBank y GFF3 desde NCBI según un listado TSV."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Ejemplo: python main.py --email usuario@institucion.org "
            "-i genomas.tsv -o descargas_ncbi"
        ),
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Correo para NCBI Entrez (obligatorio por política de NCBI).",
    )
    parser.add_argument(
        "-i",
        "--input",
        default="genomas.tsv",
        type=Path,
        help="Archivo TSV de entrada.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="descargas_ncbi",
        type=Path,
        help="Directorio base para genbank/, gff3/ y el reporte de errores.",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.4,
        help="Segundos entre peticiones (sin API key, usar >= 0.34).",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="API key opcional de NCBI (mayor cuota de peticiones).",
    )
    parser.add_argument(
        "--col-genome-id",
        default="genome_id",
        help=(
            "Columna del ID de genoma: nombre en cabecera o número 1-based "
            "(1 = primera columna)."
        ),
    )
    parser.add_argument(
        "--col-genome-name",
        default="genome_name",
        help=(
            "Columna del nombre de genoma: nombre en cabecera o número 1-based."
        ),
    )
    parser.add_argument(
        "--col-accessions",
        default="genbank_accessions",
        help=(
            "Columna de accesiones GenBank (separadas por coma): nombre o "
            "número 1-based."
        ),
    )
    return parser.parse_args(argv)


def descargar_genbank(accession: str, outdir: str) -> tuple[bool, str]:
    """Descarga registro GenBank (gbwithparts) desde Entrez nucleotide.

    Args:
        accession: Identificador NCBI nucleotide.
        outdir: Directorio donde guardar ``{accession}.gb``.

    Returns:
        Par ``(éxito, mensaje)``. Si el archivo ya existía y no está vacío,
        el mensaje indica que se omitió la descarga.
    """
    outfile = os.path.join(outdir, f"{accession}.gb")
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        return True, "Ya existía"
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gbwithparts",
            retmode="text",
        )
        data = handle.read()
        handle.close()
        if not data or "Error" in data[:200]:
            return False, "Respuesta vacía o error"
        with open(outfile, "w", encoding="utf-8") as f:
            f.write(data)
        return True, "OK"
    except Exception as e:
        return False, f"GenBank error: {e}"


def descargar_gff3(accession: str, outdir: str) -> tuple[bool, str]:
    """Descarga anotación GFF3 vía el visor sviewer de NCBI.

    Args:
        accession: Identificador NCBI nucleotide (nuccore).
        outdir: Directorio donde guardar ``{accession}.gff3``.

    Returns:
        Par ``(éxito, mensaje)``.
    """
    outfile = os.path.join(outdir, f"{accession}.gff3")
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        return True, "Ya existía"
    url = (
        "https://eutils.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"
        f"?db=nuccore&report=gff3&id={accession}"
    )
    try:
        r = requests.get(url, timeout=60)
        if r.status_code != 200 or len(r.text) < 50:
            return False, f"HTTP {r.status_code} o respuesta vacía"
        if r.text.lstrip().startswith("<"):
            return False, "Respuesta HTML (posible ID inexistente)"
        with open(outfile, "w", encoding="utf-8") as f:
            f.write(r.text)
        return True, "OK"
    except Exception as e:
        return False, f"GFF3 error: {e}"


def verificar_existencia(accession: str) -> bool:
    """Comprueba si hay al menos un resultado en nucleotide para el término.

    Args:
        accession: Término de búsqueda (normalmente el mismo accession).

    Returns:
        ``True`` si ``esearch`` devuelve ``Count`` > 0.
    """
    try:
        handle = Entrez.esearch(db="nucleotide", term=accession)
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"]) > 0
    except Exception:
        return False


def main(argv: list[str] | None = None) -> int:
    """Configura Entrez, recorre el TSV y descarga GenBank y GFF3 por accession.

    Args:
        argv: Argumentos de línea de comandos (por defecto ``sys.argv[1:]``).

    Returns:
        Código de salida 0.
    """
    args = parse_args(argv)
    Entrez.email = args.email
    Entrez.api_key = args.api_key

    output_dir = os.fspath(args.output_dir)
    genbank_dir = os.path.join(output_dir, "genbank")
    gff_dir = os.path.join(output_dir, "gff3")
    report_file = os.path.join(output_dir, "reporte_errores.tsv")

    Path(genbank_dir).mkdir(parents=True, exist_ok=True)
    Path(gff_dir).mkdir(parents=True, exist_ok=True)

    input_path = args.input
    if not input_path.is_file():
        print(f"Error: no existe el archivo de entrada: {input_path}", file=sys.stderr)
        return 1

    df = pd.read_csv(input_path, sep="\t")
    try:
        col_genome_id = resolve_column_label(args.col_genome_id, df.columns)
        col_genome_name = resolve_column_label(args.col_genome_name, df.columns)
        col_accessions = resolve_column_label(args.col_accessions, df.columns)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    delay = args.delay
    errores: list[list[str]] = []

    for _, row in df.iterrows():
        genome_id = row[col_genome_id]
        genome_name = row[col_genome_name]
        accessions = str(row[col_accessions]).split(",")
        accessions = [a.strip() for a in accessions if a.strip()]

        print(f"\n>>> {genome_id} ({genome_name})")

        for acc in accessions:
            print(f"   - {acc} ...", end=" ")

            if not verificar_existencia(acc):
                msg = "Accession NO encontrado en NCBI"
                print(f"❌ {msg}")
                errores.append([genome_id, acc, "existencia", msg])
                time.sleep(delay)
                continue

            ok_gb, msg_gb = descargar_genbank(acc, genbank_dir)
            time.sleep(delay)

            ok_gff, msg_gff = descargar_gff3(acc, gff_dir)
            time.sleep(delay)

            if ok_gb and ok_gff:
                print("✅ GenBank + GFF3")
            else:
                print(f"⚠️  GB:{msg_gb} | GFF:{msg_gff}")
                if not ok_gb:
                    errores.append([genome_id, acc, "genbank", msg_gb])
                if not ok_gff:
                    errores.append([genome_id, acc, "gff3", msg_gff])

    if errores:
        err_df = pd.DataFrame(
            errores,
            columns=["genome_id", "accession", "formato", "error"],
        )
        err_df.to_csv(report_file, sep="\t", index=False)
        print(f"\n⚠️  Reporte de errores guardado en: {report_file}")
    else:
        print("\n✅ Todas las descargas completadas sin errores.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

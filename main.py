"""Descarga GenBank, GFF3 y FASTA desde NCBI para genomas listados en un TSV.

Lee un archivo tabular con identificadores de genoma, nombres y accesiones
GenBank (separadas por comas). Para cada accession verifica existencia en
NCBI, descarga ``.gb.gz``, ``.gff3.gz`` y ``.fna.gz`` y genera reportes de
errores e integridad.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
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
            "Descarga archivos GenBank, GFF3 y FASTA desde NCBI según un TSV."
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
        "--verify-checksum",
        action="store_true",
        help="Activa verificación de checksum remoto cuando esté disponible.",
    )
    parser.add_argument(
        "--checksum-source",
        choices=("auto", "headers", "none"),
        default="auto",
        help="Fuente para checksum remoto. 'auto' intenta headers HTTP.",
    )
    parser.add_argument(
        "--checksum-timeout",
        type=float,
        default=30.0,
        help="Timeout para intentos de checksum remoto en segundos.",
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


def write_gzip_text(path: str, text: str) -> None:
    """Escribe texto comprimido en gzip con encabezado determinista."""
    data = text.encode("utf-8")
    compressed = gzip.compress(data, mtime=0)
    with open(path, "wb") as f:
        f.write(compressed)


def get_remote_checksum_from_headers(headers: requests.structures.CaseInsensitiveDict) -> tuple[str | None, str]:
    """Extrae checksum remoto de encabezados HTTP cuando sea posible.

    Returns:
        Par ``(checksum, tipo)``. Tipo en ``md5`` o ``sha256``.
    """
    content_md5 = headers.get("Content-MD5")
    if content_md5:
        return content_md5.strip(), "md5"

    etag = headers.get("ETag")
    if etag:
        candidate = etag.strip().strip('"')
        if re.fullmatch(r"[0-9a-fA-F]{32}", candidate):
            return candidate.lower(), "md5"
        if re.fullmatch(r"[0-9a-fA-F]{64}", candidate):
            return candidate.lower(), "sha256"

    return None, "none"


def checksum_file(path: str, algorithm: str = "sha256") -> str:
    """Calcula checksum para un archivo."""
    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def evaluar_checksum(
    file_path: str,
    checksum_enabled: bool,
    remote_checksum: str | None,
    remote_type: str,
) -> tuple[str, str, str]:
    """Evalúa checksum remoto y local para reportes.

    Args:
        file_path: Ruta de archivo descargado.
        checksum_enabled: Si se pidió verificación remota.
        remote_checksum: Valor remoto cuando existe.
        remote_type: Tipo de checksum remoto (md5, sha256, none).

    Returns:
        ``(status, local_sha256, remote)``.
    """
    local_sha256 = checksum_file(file_path, "sha256")
    remote_str = remote_checksum or ""
    if not checksum_enabled:
        return "local-only", local_sha256, remote_str
    if not remote_checksum or remote_type == "none":
        return "unavailable", local_sha256, remote_str
    if remote_type == "sha256":
        return ("match" if local_sha256 == remote_checksum.lower() else "mismatch"), local_sha256, remote_str
    local_md5 = checksum_file(file_path, "md5")
    return ("match" if local_md5 == remote_checksum.lower() else "mismatch"), local_sha256, remote_str


def descargar_genbank(
    accession: str,
    outdir: str,
    checksum_enabled: bool,
) -> tuple[bool, str, str, str]:
    """Descarga registro GenBank (gbwithparts) y lo guarda en ``.gb.gz``."""
    outfile = os.path.join(outdir, f"{accession}.gb.gz")
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        status, local_sha, remote = evaluar_checksum(outfile, checksum_enabled, None, "none")
        return True, "Ya existía", status, f"{local_sha}|{remote}"
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
            return False, "Respuesta vacía o error", "unavailable", ""
        write_gzip_text(outfile, data)
        status, local_sha, remote = evaluar_checksum(outfile, checksum_enabled, None, "none")
        return True, "OK", status, f"{local_sha}|{remote}"
    except Exception as e:
        return False, f"GenBank error: {e}", "unavailable", ""


def descargar_gff3(
    accession: str,
    outdir: str,
    checksum_enabled: bool,
    checksum_source: str,
    checksum_timeout: float,
) -> tuple[bool, str, str, str]:
    """Descarga anotación GFF3 vía el visor sviewer de NCBI en ``.gff3.gz``."""
    del checksum_timeout  # Reservado para futuras fuentes remotas distintas.
    outfile = os.path.join(outdir, f"{accession}.gff3.gz")
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        status, local_sha, remote = evaluar_checksum(outfile, checksum_enabled, None, "none")
        return True, "Ya existía", status, f"{local_sha}|{remote}"
    url = (
        "https://eutils.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"
        f"?db=nuccore&report=gff3&id={accession}"
    )
    try:
        r = requests.get(url, timeout=60)
        if r.status_code != 200 or len(r.text) < 50:
            return False, f"HTTP {r.status_code} o respuesta vacía", "unavailable", ""
        if r.text.lstrip().startswith("<"):
            return False, "Respuesta HTML (posible ID inexistente)", "unavailable", ""
        write_gzip_text(outfile, r.text)
        remote_checksum = None
        remote_type = "none"
        if checksum_enabled and checksum_source in {"auto", "headers"}:
            remote_checksum, remote_type = get_remote_checksum_from_headers(r.headers)
        status, local_sha, remote = evaluar_checksum(
            outfile, checksum_enabled, remote_checksum, remote_type
        )
        return True, "OK", status, f"{local_sha}|{remote}"
    except Exception as e:
        return False, f"GFF3 error: {e}", "unavailable", ""


def descargar_fasta(
    accession: str,
    outdir: str,
    checksum_enabled: bool,
) -> tuple[bool, str, str, str]:
    """Descarga secuencia FASTA desde Entrez y la guarda en ``.fna.gz``."""
    outfile = os.path.join(outdir, f"{accession}.fna.gz")
    if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        status, local_sha, remote = evaluar_checksum(outfile, checksum_enabled, None, "none")
        return True, "Ya existía", status, f"{local_sha}|{remote}"
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="fasta",
            retmode="text",
        )
        data = handle.read()
        handle.close()
        if not data or not data.lstrip().startswith(">"):
            return False, "Respuesta vacía o no FASTA", "unavailable", ""
        write_gzip_text(outfile, data)
        status, local_sha, remote = evaluar_checksum(outfile, checksum_enabled, None, "none")
        return True, "OK", status, f"{local_sha}|{remote}"
    except Exception as e:
        return False, f"FASTA error: {e}", "unavailable", ""


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
    """Descarga GenBank, GFF3 y FASTA comprimidos por accession.

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
    fasta_dir = os.path.join(output_dir, "fasta")
    report_file = os.path.join(output_dir, "reporte_errores.tsv")
    checksum_report_file = os.path.join(output_dir, "reporte_checksums.tsv")

    Path(genbank_dir).mkdir(parents=True, exist_ok=True)
    Path(gff_dir).mkdir(parents=True, exist_ok=True)
    Path(fasta_dir).mkdir(parents=True, exist_ok=True)

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
    checksum_rows: list[list[str]] = []

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

            ok_gb, msg_gb, gb_status, gb_checks = descargar_genbank(
                acc,
                genbank_dir,
                args.verify_checksum,
            )
            time.sleep(delay)

            ok_gff, msg_gff, gff_status, gff_checks = descargar_gff3(
                acc,
                gff_dir,
                args.verify_checksum,
                args.checksum_source,
                args.checksum_timeout,
            )
            time.sleep(delay)

            ok_fna, msg_fna, fna_status, fna_checks = descargar_fasta(
                acc,
                fasta_dir,
                args.verify_checksum,
            )
            time.sleep(delay)

            for formato, status, checks in (
                ("genbank", gb_status, gb_checks),
                ("gff3", gff_status, gff_checks),
                ("fasta", fna_status, fna_checks),
            ):
                local_sha, remote = (checks.split("|", 1) if "|" in checks else ("", ""))
                checksum_rows.append(
                    [genome_id, acc, formato, status, local_sha, remote]
                )

            if ok_gb and ok_gff and ok_fna:
                print(
                    "✅ GenBank + GFF3 + FASTA "
                    f"[chk: gb={gb_status}, gff={gff_status}, fa={fna_status}]"
                )
            else:
                print(
                    "⚠️  "
                    f"GB:{msg_gb} | GFF:{msg_gff} | FASTA:{msg_fna} "
                    f"[chk: gb={gb_status}, gff={gff_status}, fa={fna_status}]"
                )
                if not ok_gb:
                    errores.append([genome_id, acc, "genbank", msg_gb])
                if not ok_gff:
                    errores.append([genome_id, acc, "gff3", msg_gff])
                if not ok_fna:
                    errores.append([genome_id, acc, "fasta", msg_fna])

    if errores:
        err_df = pd.DataFrame(
            errores,
            columns=["genome_id", "accession", "formato", "error"],
        )
        err_df.to_csv(report_file, sep="\t", index=False)
        print(f"\n⚠️  Reporte de errores guardado en: {report_file}")
    else:
        print("\n✅ Todas las descargas completadas sin errores.")

    if checksum_rows:
        checksum_df = pd.DataFrame(
            checksum_rows,
            columns=[
                "genome_id",
                "accession",
                "formato",
                "checksum_status",
                "local_sha256",
                "remote_checksum",
            ],
        )
        checksum_df.to_csv(checksum_report_file, sep="\t", index=False)
        print(f"ℹ️  Reporte de checksums guardado en: {checksum_report_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Check AlphaFold DB (AFDB) coverage for proteins listed in a FASTA file.

What this script does:
1. Reads a FASTA file containing protein sequences
2. Extracts UniProt accessions from FASTA headers
3. Checks whether each accession has a structure in AlphaFold DB
4. Writes a TSV report with the results

Why this version:
- Fully compatible with Python 3.12
- Uses ONLY the Python standard library
- No external dependencies (no requests, no tqdm)
"""

# --------------------------------------------------
# Imports (standard library only)
# --------------------------------------------------

from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
import argparse
import time


# --------------------------------------------------
# AlphaFold DB URL template
# --------------------------------------------------
# AFDB stores one main structure per UniProt accession:
# AF-{UniProtAC}-F1-model_v4.cif
#
# If this URL exists → the protein is in AlphaFold DB
AFDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{ac}-F1-model_v6.cif"



# --------------------------------------------------
# FASTA header parsing
# --------------------------------------------------
def extract_uniprot_ac(header: str) -> str:
    """
    Extract a UniProt accession from a FASTA header line.

    Supported header formats:
    >Q9AD92
    >Q9AD92 some description
    >sp|Q9AD92|CHPH_STRCO
    >tr|A0A2K6V5L6|...

    Strategy:
    - Remove '>' character
    - If the header contains '|', take the second field
    - Otherwise, take the first whitespace-separated token
    """
    header = header.lstrip(">").strip()

    if "|" in header:
        parts = header.split("|")
        if len(parts) >= 2:
            return parts[1]

    return header.split()[0]


# --------------------------------------------------
# Read all accessions from FASTA
# --------------------------------------------------
def read_accessions(fasta: Path) -> list[str]:
    """
    Read a FASTA file and return a sorted list of UNIQUE
    UniProt accessions found in the headers.
    """
    accs = []

    with fasta.open() as f:
        for line in f:
            if line.startswith(">"):
                accs.append(extract_uniprot_ac(line))

    # Remove duplicates and sort
    return sorted(set(accs))


# --------------------------------------------------
# Check if an accession exists in AlphaFold DB
# --------------------------------------------------
def afdb_exists(accession: str) -> bool:
    """
    Check whether a UniProt accession exists in AlphaFold DB.

    How this works:
    - Send an HTTP HEAD request to the expected AFDB URL
    - HEAD is fast and downloads no file content
    - If the server responds → structure exists
    """
    url = AFDB_URL.format(ac=accession)
    req = Request(url, method="HEAD")

    try:
        with urlopen(req, timeout=10):
            return True
    except HTTPError as e:
        # Some servers return 200 as an exception
        return e.code == 200
    except URLError:
        # Connection problems or missing file
        return False


# --------------------------------------------------
# Main program logic
# --------------------------------------------------
def main():
    """
    Main entry point of the script.
    Handles:
    - Command-line arguments
    - Running AFDB checks
    - Writing output
    """

    # -----------------------------
    # Command-line arguments
    # -----------------------------
    parser = argparse.ArgumentParser(
        description="Check AlphaFold DB coverage for a FASTA file"
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Input FASTA file with protein sequences",
    )
    parser.add_argument(
        "--out",
        default=Path("afdb_coverage.tsv"),
        type=Path,
        help="Output TSV file (default: afdb_coverage.tsv)",
    )

    args = parser.parse_args()

    # -----------------------------
    # Read UniProt accessions
    # -----------------------------
    accessions = read_accessions(args.fasta)
    print(f"Found {len(accessions)} unique UniProt accessions")

    # -----------------------------
    # Open output file
    # -----------------------------
    with args.out.open("w") as out:
        # Write header
        out.write("uniprot_ac\tpresent_in_afdb\tmodel_url\n")

        # -----------------------------
        # Check each accession
        # -----------------------------
        for i, ac in enumerate(accessions, start=1):
            present = afdb_exists(ac)
            url = AFDB_URL.format(ac=ac)

            out.write(f"{ac}\t{'yes' if present else 'no'}\t{url}\n")

            # Small delay to avoid stressing AFDB servers
            time.sleep(0.1)

            # Progress update every 25 proteins
            if i % 25 == 0:
                print(f"Checked {i}/{len(accessions)}")

    print(f"\nAFDB coverage report written to: {args.out}")


# --------------------------------------------------
# Script entry point
# --------------------------------------------------
if __name__ == "__main__":
    main()

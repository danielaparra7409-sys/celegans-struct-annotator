#!/usr/bin/env python3
"""
Script local para:
1) Leer entradas desde FASTA
2) Llamar al endpoint de AlphaFold2 de NVIDIA
3) Guardar JSON y PDB por Entry en predicciones/
"""

import os
import json
import time
import requests
from pathlib import Path
from typing import Dict, Optional

print(">>> estructuras.py se está ejecutando <<<")

# --------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------
FASTA_SHORT = "microdomains_short.fasta"
FASTA_LONG = "microdomains_long.fasta"

OUT_DIR = Path("predicciones")
OUT_DIR.mkdir(exist_ok=True)

ALPHAFOLD_URL = "https://health.api.nvidia.com/v1/alphafold2"
STATUS_URL    = "https://health.api.nvidia.com/v1/status"


# --------------------------------------------------------------------
# JSON → PDB
# --------------------------------------------------------------------
def json_to_pdb(json_file: str, pdb_file: str) -> None:
    with open(json_file, "r") as f:
        json_data = json.load(f)

    if isinstance(json_data, list):
        with open(pdb_file, "w") as f:
            for line in json_data:
                f.write(line)
    else:
        raise ValueError("Formato JSON inesperado.")


# --------------------------------------------------------------------
# Leer FASTA como {entry: secuencia}
# --------------------------------------------------------------------
def load_fasta_sequences(path: str) -> Dict[str, str]:
    seqs = {}
    current_id = None
    current_seq = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                # Guardar secuencia previa
                if current_id and current_seq:
                    seqs[current_id] = "".join(current_seq)

                # Obtener ID limpio
                header = line[1:].strip()
                current_id = header.split()[0]
                current_seq = []
            else:
                current_seq.append("".join([c for c in line if c.isalpha()]))

        # Guardar último registro
        if current_id and current_seq:
            seqs[current_id] = "".join(current_seq)

    return seqs


# --------------------------------------------------------------------
# Llamado a AlphaFold2
# --------------------------------------------------------------------
def run_alphafold_for_entry(entry: str, sequence: str, key: str) -> Optional[Path]:
    json_path = OUT_DIR / f"{entry}.json"

    headers = {
        "content-type": "application/json",
        "Authorization": f"Bearer {key}",
        "NVCF-POLL-SECONDS": "300",
    }

    data = {
        "sequence": sequence,
        "algorithm": "mmseqs2",
        "e_value": 0.0001,
        "iterations": 1,
        "databases": ["small_bfd"],
        "relax_prediction": False,
        "skip_template_search": True,
    }

    print(f"\n[INFO] Enviando petición para Entry: {entry} ({len(sequence)} aa)")

    try:
        response = requests.post(ALPHAFOLD_URL, headers=headers, json=data)
    except Exception as e:
        print(f"[ERROR] Conexión fallida: {e}")
        return None

    # --- Respuesta directa
    if response.status_code == 200:
        json_path.write_text(response.text)
        print(f"[OK] Respuesta guardada en {json_path}")
        return json_path

    # --- Caso 202 → Polling
    if response.status_code == 202:
        print("[INFO] Esperando resultado (202)…")
        req_id = response.headers.get("nvcf-reqid")
        if not req_id:
            print("[ERROR] No se recibió nvcf-reqid.")
            return None

        while True:
            status = requests.get(f"{STATUS_URL}/{req_id}", headers=headers)

            if status.status_code == 202:
                print("[INFO] Procesando… reintentando en 8s")
                time.sleep(8)
                continue

            json_path.write_text(status.text)
            print(f"[OK] Respuesta final guardada en {json_path}")
            return json_path

    print(f"[ERROR] HTTP {response.status_code}: {response.text}")
    return None


# --------------------------------------------------------------------
# MODO INTERACTIVO
# --------------------------------------------------------------------
def main_interactivo():
    print("\n=== MODO INTERACTIVO ===")
    print("• Usa FASTA short y long")
    print("• Ingresas un Entry")
    print("• Produce JSON + PDB")
    print("• Enter vacío → salir\n")

    key = os.getenv("NVCF_RUN_KEY") or input("Pega tu Run Key de NVIDIA: ").strip()
    if not key:
        print("ERROR: No se proporcionó API Key.")
        return

    fasta_short = load_fasta_sequences(FASTA_SHORT)
    fasta_long = load_fasta_sequences(FASTA_LONG)

    while True:
        print("\n----------------------------------------")
        entry = input("Entry (Enter vacío para salir): ").strip()

        if entry == "":
            print("Fin del modo interactivo.")
            break

        seq = fasta_short.get(entry) or fasta_long.get(entry)

        if seq is None:
            print(f"[AVISO] No encontré secuencia para '{entry}'.")
            continue

        print(f"[INFO] Secuencia encontrada: {len(seq)} aa")

        if len(seq) < 50:
            print(f"[AVISO] {entry} tiene {len(seq)} aa → NVIDIA no acepta < 50 aa.")
            continue

        json_path = run_alphafold_for_entry(entry, seq, key)

        if json_path is None:
            print(f"[ERROR] Falló {entry}")
            continue

        pdb_path = OUT_DIR / f"{entry}.pdb"
        try:
            json_to_pdb(str(json_path), str(pdb_path))
            print(f"[OK] PDB generado: {pdb_path}")
        except Exception as e:
            print(f"[ERROR] Fallo JSON→PDB: {e}")


# --------------------------------------------------------------------
if __name__ == "__main__":
    main_interactivo()

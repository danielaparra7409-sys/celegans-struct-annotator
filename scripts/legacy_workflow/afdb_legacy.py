#!/usr/bin/env python3
"""
afdb_domains_pipeline_v1.py

Pipeline reproducible:
1) Descarga modelos AlphaFold DB por UniProt accession (AFDB) -> mmCIF
2) (Opcional) Convierte mmCIF -> PDB
3) Recorta dominios (rango de residuos) desde mmCIF/PDB completos -> archivos por dominio
4) Genera links Mol* para visualizar dominios en 3D (sin subir manualmente)

Requisitos:
    pip install aiohttp pandas gemmi tqdm

Uso rápido:
    python afdb_domains_pipeline_v1.py download --fasta proteins_full.fasta
    python afdb_domains_pipeline_v1.py slice --domains domains.tsv
    python afdb_domains_pipeline_v1.py links --domains domains.tsv --base-url http://localhost:8000
"""

from __future__ import annotations

import argparse
import asyncio
import re
from pathlib import Path
from urllib.parse import quote

import aiohttp
import pandas as pd
import gemmi
from tqdm import tqdm


AFDB_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{ac}-F1-model_v4.cif"


# ---------- Utils ----------
def parse_uniprot_accessions_from_fasta(fasta_path: Path) -> list[str]:
    accessions: list[str] = []
    with fasta_path.open() as f:
        for line in f:
            if line.startswith(">"):
                # toma el primer token del header como ID
                token = line[1:].strip().split()[0]
                token = token.split("|")[0]  # por si viene "sp|AC|NAME"
                # filtra uniprot-like (simple)
                if re.match(r"^[A-Z0-9]{6,10}$", token):
                    accessions.append(token)
                else:
                    accessions.append(token)  # igual lo intentamos; si no existe dará 404
    # únicos preservando orden
    seen = set()
    uniq = []
    for a in accessions:
        if a not in seen:
            uniq.append(a)
            seen.add(a)
    return uniq


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


async def download_one(session: aiohttp.ClientSession, url: str, out_path: Path, timeout_s: int = 120) -> tuple[bool, str]:
    try:
        async with session.get(url, timeout=aiohttp.ClientTimeout(total=timeout_s)) as r:
            if r.status == 200:
                data = await r.read()
                out_path.write_bytes(data)
                return True, "ok"
            elif r.status == 404:
                return False, "404_not_found"
            else:
                text = await r.text()
                return False, f"http_{r.status}:{text[:120]}"
    except Exception as e:
        return False, f"error:{type(e).__name__}:{e}"


async def download_models(accessions: list[str], out_dir: Path, concurrency: int = 25) -> pd.DataFrame:
    ensure_dir(out_dir)
    sem = asyncio.Semaphore(concurrency)

    results = []

    async with aiohttp.ClientSession() as session:
        async def task(ac: str):
            async with sem:
                out_path = out_dir / f"{ac}.cif"
                if out_path.exists() and out_path.stat().st_size > 1000:
                    results.append((ac, True, "cached", str(out_path)))
                    return

                url = AFDB_URL_TEMPLATE.format(ac=ac)
                ok, msg = await download_one(session, url, out_path)
                if ok:
                    results.append((ac, True, msg, str(out_path)))
                else:
                    # si falló, borramos archivo vacío si existe
                    if out_path.exists() and out_path.stat().st_size < 1000:
                        out_path.unlink(missing_ok=True)
                    results.append((ac, False, msg, str(out_path)))

        tasks = [task(ac) for ac in accessions]

        for fut in tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="Downloading AFDB models"):
            await fut

    df = pd.DataFrame(results, columns=["uniprot_ac", "success", "status", "path"])
    return df


def mmcif_to_pdb(cif_path: Path, pdb_path: Path) -> None:
    st = gemmi.read_structure(str(cif_path))
    st.remove_empty_chains()
    # escribe PDB
    pdb_path.write_text(st.make_pdb_string())


def slice_domain_from_cif(full_cif: Path, out_cif: Path, start: int, end: int, chain_id: str | None = None) -> None:
    """
    Recorta por rango de residuo (seqid.num). Si chain_id es None, recorta todas las cadenas.
    """
    st = gemmi.read_structure(str(full_cif))
    st.remove_empty_chains()

    for model in st:
        for chain in list(model):  # copia para poder editar
            if chain_id and chain.name != chain_id:
                # si se especifica cadena, borra otras
                model.remove_chain(chain.name)
                continue

            kept = []
            for res in chain:
                # seqid.num puede ser None en casos raros
                num = res.seqid.num
                if num is None:
                    continue
                if start <= int(num) <= end:
                    kept.append(res)

            # reemplaza residuos
            chain.clear()
            for r in kept:
                chain.add_residue(r)

            # si quedó vacía, eliminar cadena
            if len(chain) == 0:
                model.remove_chain(chain.name)

    # escribir mmCIF
    doc = st.make_mmcif_document()
    doc.write_file(str(out_cif))


def slice_domains(domains_tsv: Path, models_dir: Path, out_dir: Path, also_pdb: bool = True, chain_col: str | None = None) -> pd.DataFrame:
    ensure_dir(out_dir)
    df = pd.read_csv(domains_tsv, sep="\t")

    required = {"uniprot_ac", "domain_id", "start", "end"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"domains.tsv debe tener columnas: {sorted(required)}. Faltan: {sorted(missing)}")

    outputs = []
    for row in tqdm(df.itertuples(index=False), total=len(df), desc="Slicing domains"):
        ac = str(getattr(row, "uniprot_ac"))
        dom = str(getattr(row, "domain_id"))
        start = int(getattr(row, "start"))
        end = int(getattr(row, "end"))
        chain_id = None
        if chain_col and hasattr(row, chain_col):
            val = getattr(row, chain_col)
            chain_id = None if pd.isna(val) else str(val)

        full_cif = models_dir / f"{ac}.cif"
        if not full_cif.exists():
            outputs.append((dom, ac, False, "missing_full_model", "", ""))
            continue

        out_cif = out_dir / f"{dom}.cif"
        out_pdb = out_dir / f"{dom}.pdb"

        try:
            slice_domain_from_cif(full_cif, out_cif, start, end, chain_id=chain_id)
            if also_pdb:
                mmcif_to_pdb(out_cif, out_pdb)
            outputs.append((dom, ac, True, "ok", str(out_cif), str(out_pdb if also_pdb else "")))
        except Exception as e:
            outputs.append((dom, ac, False, f"error:{type(e).__name__}:{e}", "", ""))

    out_df = pd.DataFrame(outputs, columns=["domain_id", "uniprot_ac", "success", "status", "cif_path", "pdb_path"])
    return out_df


def make_links(domains_tsv: Path, out_links: Path, base_url: str, fmt: str = "pdb") -> None:
    """
    Genera links Mol* por dominio.
    base_url: URL donde estás sirviendo los archivos (ej. http://localhost:8000)
    fmt: "pdb" o "mmcif"
    """
    df = pd.read_csv(domains_tsv, sep="\t")
    fmt = "pdb" if fmt.lower() in ("pdb",) else "mmcif"
    ext = "pdb" if fmt == "pdb" else "cif"
    molfmt = "pdb" if fmt == "pdb" else "mmcif"

    def link_for(dom_id: str) -> str:
        file_url = f"{base_url.rstrip('/')}/{dom_id}.{ext}"
        return f"https://molstar.org/viewer/?structure-url={quote(file_url)}&structure-url-format={molfmt}"

    df["molstar_url"] = df["domain_id"].astype(str).apply(link_for)
    df[["uniprot_ac", "domain_id", "start", "end", "molstar_url"]].to_csv(out_links, sep="\t", index=False)


# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_dl = sub.add_parser("download", help="Descarga modelos AFDB por UniProt AC desde un FASTA")
    ap_dl.add_argument("--fasta", type=Path, required=True)
    ap_dl.add_argument("--out", type=Path, default=Path("models_cif"))
    ap_dl.add_argument("--concurrency", type=int, default=25)
    ap_dl.add_argument("--report", type=Path, default=Path("download_report.tsv"))

    ap_slice = sub.add_parser("slice", help="Recorta dominios desde modelos completos (mmCIF) usando domains.tsv")
    ap_slice.add_argument("--domains", type=Path, required=True)
    ap_slice.add_argument("--models", type=Path, default=Path("models_cif"))
    ap_slice.add_argument("--out", type=Path, default=Path("domains_out"))
    ap_slice.add_argument("--no-pdb", action="store_true", help="No generar PDB (solo CIF)")
    ap_slice.add_argument("--report", type=Path, default=Path("slice_report.tsv"))
    ap_slice.add_argument("--chain-col", type=str, default=None, help="Nombre de columna opcional con chain_id")

    ap_links = sub.add_parser("links", help="Genera TSV de links Mol* para visualizar dominios")
    ap_links.add_argument("--domains", type=Path, required=True)
    ap_links.add_argument("--out", type=Path, default=Path("links_domains.tsv"))
    ap_links.add_argument("--base-url", type=str, required=True, help="Ej: http://localhost:8000 (donde sirves domains_out)")
    ap_links.add_argument("--format", type=str, default="pdb", choices=["pdb", "mmcif"])

    args = ap.parse_args()

    if args.cmd == "download":
        accessions = parse_uniprot_accessions_from_fasta(args.fasta)
        df = asyncio.run(download_models(accessions, args.out, concurrency=args.concurrency))
        df.to_csv(args.report, sep="\t", index=False)
        ok = int(df["success"].sum())
        total = len(df)
        print(f"\nDescargados OK: {ok}/{total}")
        print(f"Reporte: {args.report}")

    elif args.cmd == "slice":
        df = slice_domains(
            domains_tsv=args.domains,
            models_dir=args.models,
            out_dir=args.out,
            also_pdb=(not args.no_pdb),
            chain_col=args.chain_col
        )
        df.to_csv(args.report, sep="\t", index=False)
        ok = int(df["success"].sum())
        total = len(df)
        print(f"\nDominios recortados OK: {ok}/{total}")
        print(f"Reporte: {args.report}")
        print(f"Salida dominios: {args.out}")

    elif args.cmd == "links":
        make_links(args.domains, args.out, base_url=args.base_url, fmt=args.format)
        print(f"Links generados: {args.out}")

if __name__ == "__main__":
    main()

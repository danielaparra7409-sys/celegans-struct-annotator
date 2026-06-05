from pathlib import Path
import pandas as pd
import urllib.request
import urllib.error
import json
import time

BASE = Path(
    "/mnt/c/Users/Daniela/OneDrive - Universidad Nacional de Colombia"
    "/Escritorio/Tareas/Trabajo C.elegans/CeMbio"
)

AF_API = "https://alphafold.ebi.ac.uk/api/prediction"
USER_AGENT = "celegans-struct-annotator/1.0 (research)"
DELAY = 1.5  # seconds between requests

in_table = BASE / "tables" / "nls50_swissprot_accessions_for_alphafold.tsv"
out_note = BASE / "notes" / "nls50_alphafold_pilot_fetch_summary.txt"
struct_dir = BASE / "structures" / "alphafold"
struct_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(in_table, sep="\t", dtype=str)

# Work only on priority 1+2 (strong + moderate) rows still pending — idempotent
mask = (
    df["alphafold_lookup_priority"].fillna("99").astype(int) <= 2
) & (df["alphafold_status"] == "pending")
targets = df[mask]

print(f"AlphaFold pilot fetch — priority targets (strong+moderate): {len(targets)}")
print(f"Structures will be saved to: {struct_dir}")
print()


def _get(url: str, binary: bool = False):
    req = urllib.request.Request(url)
    req.add_header("User-Agent", USER_AGENT)
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read() if binary else json.loads(resp.read())


def fetch_accession(acc: str, cif_path: Path, pdb_path: Path) -> dict:
    api_data = _get(f"{AF_API}/{acc}")
    if not api_data:
        return {"alphafold_status": "not_found", "alphafold_lookup_note": "empty API response"}

    entry = api_data[0]
    cif_url = entry["cifUrl"]
    pdb_url = entry["pdbUrl"]
    version = entry.get("latestVersion", "?")
    plddt = entry.get("globalMetricValue", "")

    time.sleep(DELAY)
    cif_path.write_bytes(_get(cif_url, binary=True))
    time.sleep(DELAY)
    pdb_path.write_bytes(_get(pdb_url, binary=True))

    note = f"AlphaFold DB v{version}, pLDDT={plddt}"
    return {
        "alphafold_status": "found",
        "alphafold_cif_url": cif_url,
        "alphafold_pdb_url": pdb_url,
        "alphafold_cif_path": str(cif_path),
        "alphafold_pdb_path": str(pdb_path),
        "alphafold_lookup_note": note,
    }


updates = {}

for _, row in targets.iterrows():
    acc = row["swissprot_accession"]
    label = row.get("swissprot_entry_from_fasta", acc) or acc
    conf = row["manual_review_confidence"]
    cif_path = struct_dir / f"AF-{acc}-F1-model_latest.cif"
    pdb_path = struct_dir / f"AF-{acc}-F1-model_latest.pdb"

    print(f"  {acc}  {label:<20}  ({conf}) ...", end=" ", flush=True)

    try:
        result = fetch_accession(acc, cif_path, pdb_path)
        updates[acc] = result
        if result["alphafold_status"] == "found":
            print(f"found  [{result['alphafold_lookup_note']}]")
        else:
            print(result["alphafold_status"])
        time.sleep(DELAY)

    except urllib.error.HTTPError as e:
        if e.code == 404:
            status_label = "not_found"
            note = "not in AlphaFold DB (404)"
        else:
            status_label = "error"
            note = f"HTTP {e.code}"
        updates[acc] = {
            "alphafold_status": status_label,
            "alphafold_cif_url": "",
            "alphafold_pdb_url": "",
            "alphafold_cif_path": "",
            "alphafold_pdb_path": "",
            "alphafold_lookup_note": note,
        }
        print(status_label)
        time.sleep(DELAY)

    except Exception as e:
        updates[acc] = {
            "alphafold_status": "error",
            "alphafold_cif_url": "",
            "alphafold_pdb_url": "",
            "alphafold_cif_path": "",
            "alphafold_pdb_path": "",
            "alphafold_lookup_note": str(e)[:200],
        }
        print(f"error: {e}")
        time.sleep(DELAY)

# Write updates back to df
for acc, upd in updates.items():
    idx = df.index[df["swissprot_accession"] == acc]
    for col, val in upd.items():
        df.loc[idx, col] = val

df.to_csv(in_table, sep="\t", index=False)

# Summary counts
found = [acc for acc, u in updates.items() if u["alphafold_status"] == "found"]
not_found = [acc for acc, u in updates.items() if u["alphafold_status"] == "not_found"]
errors = [acc for acc, u in updates.items() if u["alphafold_status"] == "error"]

with out_note.open("w", encoding="utf-8") as f:
    f.write("NLS 50 percent AlphaFold pilot fetch summary\n")
    f.write("============================================\n\n")

    f.write("Purpose\n")
    f.write("-------\n")
    f.write(
        "Retrieve AlphaFold DB structures for the 9 priority SwissProt accessions "
        "(strong + moderate confidence hits) from the NLS 50% representative pilot search.\n\n"
    )

    f.write("Parameters\n")
    f.write("----------\n")
    f.write(f"AlphaFold API: {AF_API}\n")
    f.write("Priority filter: alphafold_lookup_priority <= 2 (strong + moderate)\n")
    f.write(f"Structure output directory: {struct_dir}\n\n")

    f.write("Results\n")
    f.write("-------\n")
    f.write(f"Targets attempted: {len(updates)}\n")
    f.write(f"  Found (downloaded): {len(found)}\n")
    f.write(f"  Not in AlphaFold DB: {len(not_found)}\n")
    f.write(f"  Errors: {len(errors)}\n\n")

    if found:
        f.write("Found accessions\n")
        f.write("----------------\n")
        for acc in found:
            row = df[df["swissprot_accession"] == acc].iloc[0]
            note = updates[acc].get("alphafold_lookup_note", "")
            f.write(
                f"  {acc}  {row['swissprot_entry_from_fasta']}  "
                f"{row['swissprot_protein_name']}  "
                f"({row['manual_review_confidence']})  [{note}]\n"
            )
        f.write("\n")

    if not_found:
        f.write("Not found in AlphaFold DB\n")
        f.write("-------------------------\n")
        for acc in not_found:
            row = df[df["swissprot_accession"] == acc].iloc[0]
            f.write(
                f"  {acc}  {row['swissprot_entry_from_fasta']}  "
                f"{row['swissprot_protein_name']}  ({row['manual_review_confidence']})\n"
            )
        f.write("\n")

    if errors:
        f.write("Errors\n")
        f.write("------\n")
        for acc in errors:
            f.write(f"  {acc}  {updates[acc]['alphafold_lookup_note']}\n")
        f.write("\n")

    f.write("Caution\n")
    f.write("-------\n")
    f.write(
        "Structures retrieved are for the SwissProt homologs, not the original CeMbio proteins.\n"
        "AlphaFold DB coverage is limited to proteins in UniProt; phage and low-PE proteins "
        "may be absent. Not-found accessions should be considered for ESMFold or structure "
        "prediction from the CeMbio representative sequence directly.\n"
    )

print()
print(f"Updated table : {in_table}")
print(f"Summary note  : {out_note}")
print()
print(f"Found: {len(found)}  |  Not found: {len(not_found)}  |  Errors: {len(errors)}")
if found:
    print(f"  Downloaded: {', '.join(found)}")
if not_found:
    print(f"  Missing:    {', '.join(not_found)}")

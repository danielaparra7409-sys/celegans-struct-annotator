import csv
import requests
from pathlib import Path

input_tsv = "afdb_coverage.tsv"
output_dir = Path("afdb_cif")
output_dir.mkdir(exist_ok=True)

with open(input_tsv) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["present_in_afdb"] == "yes":
            url = row["model_url"]
            filename = url.split("/")[-1]
            out_path = output_dir / filename

            print(f"Downloading {filename}...")
            r = requests.get(url)
            out_path.write_bytes(r.content)

print("Done.")

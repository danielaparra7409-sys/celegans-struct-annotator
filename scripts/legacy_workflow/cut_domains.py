import os
import re
import pandas as pd
from Bio import SeqIO

# === CONFIGURACIÓN DE RUTAS ===
excel_path = r'C:/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/BD_completed.xlsx'
fasta_files = [
    r'C:/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/Domains/TL_secreted.fasta',
    r'C:/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/Domains/PL_secreted.fasta',
    r'C:/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/Domains/Homology_secreted.fasta'
]

output_dir = r'C:/Users/Daniela/OneDrive - Universidad Nacional de Colombia/Escritorio/Tareas/Trabajo C.elegans/domains'
os.makedirs(output_dir, exist_ok=True)
output_fasta = os.path.join(output_dir, 'microdomains.fasta')
output_txt = os.path.join(output_dir, 'microdomains.txt')

# === CARGA DE DATOS ===
print("Cargando archivo Excel...")
df = pd.read_excel(excel_path)

print("Cargando archivos FASTA...")
sequences = {}
for fasta_path in fasta_files:
    for record in SeqIO.parse(fasta_path, "fasta"):
        acc = record.id.split("|")[1] if "|" in record.id else record.id  # Ajustar según formato 'sp|ID|NAME'
        sequences[acc] = str(record.seq)

print(f"Total secuencias cargadas: {len(sequences)}")

# === PROCESAMIENTO ===
results = []

for _, row in df.iterrows():
    entry = str(row['Entry']).strip()
    domains_text = str(row.get('Domain [FT]', ''))
    length = row.get('Length', '')

    if entry in sequences and 'DOMAIN' in domains_text:
        matches = re.findall(r'DOMAIN (/d+)/./.(/d+)', domains_text)
        for start, end in matches:
            start, end = int(start), int(end)
            seq = sequences[entry]
            if end <= len(seq):
                subseq = seq[start - 1:end]
                domain_len = len(subseq)
                header = f">{entry} Domain:{start}-{end} Len_Domain:{domain_len} Len_Protein:{length}"

                results.append((header, subseq))

# === GUARDADO DE RESULTADOS ===
print(f"Guardando {len(results)} microdominios...")
with open(output_fasta, 'w') as f_fasta, open(output_txt, 'w') as f_txt:
    for header, seq in results:
        f_fasta.write(header + "/n")
        for i in range(0, len(seq), 70):
            f_fasta.write(seq[i:i + 70] + "/n")
        f_txt.write(header + "/n" + seq + "/n/n")

print(f"Archivos generados:/nFASTA → {output_fasta}/nTXT → {output_txt}")
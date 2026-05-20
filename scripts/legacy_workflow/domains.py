# Paso 1: Cargar y reconocer campos de FASTA + Excel 
# Requisitos: pandas, biopython (pip install pandas biopython)

import os, re
import pandas as pd
from Bio import SeqIO

base_dir = r"C:\Users\Daniela\OneDrive - Universidad Nacional de Colombia\Escritorio\Tareas\Trabajo C.elegans\Domains"
excel_path = os.path.join(base_dir, "BD_completed.xlsx")  # está en la misma carpeta Domains


fasta_paths = {
    "TL_secreted":           os.path.join(base_dir, "TL_secreted.fasta"),
    "PL_secreted":           os.path.join(base_dir, "PL_secreted.fasta"),
    "Homology_experimental": os.path.join(base_dir, "Homology_experimental.fasta"),
    "Homology_any":          os.path.join(base_dir, "Homology_any.fasta")
}

# encabezados UniProt (sp|ACC|ID DESCRIPCION...)
_header_re = re.compile(r'^(?:\w+\|)?(?P<acc>[^|>\s]+)\|(?P<id>[^ \t|]+)\s*(?P<desc>.*)$')

def parse_fasta_to_df(fpath: str, source_label: str) -> pd.DataFrame:
    rows = []
    for rec in SeqIO.parse(fpath, "fasta"):
        raw_header = rec.description.strip()
        acc, entry_name, desc = None, None, None
        m = _header_re.match(raw_header)
        if m:
            acc        = m.group("acc")        # <- Entry (accesión UniProt)
            entry_name = m.group("id")         # <- Entry name (e.g., LYSO_STAAU)
            desc       = m.group("desc") or ""
        else:
            # Si no es formato UniProt, usa el id entero como acc
            acc, entry_name, desc = rec.id, rec.id, raw_header

        rows.append({
            "SourceFile": source_label,
            "RawHeader": raw_header,
            "Entry": acc,               # clave para enlazar con el Excel
            "EntryName": entry_name,
            "Description": desc,
            "Length": len(rec.seq),
            "Sequence": str(rec.seq)
        })
    return pd.DataFrame(rows)

# ----- LECTURA DE FASTA -----
fasta_tables = {}
for label, path in fasta_paths.items():
    if not os.path.exists(path):
        print(f"⚠️  No se encontró FASTA: {path}")
        continue
    df_fa = parse_fasta_to_df(path, label)
    fasta_tables[label] = df_fa
    print(f"\n✅ Cargado {label}: {len(df_fa)} secuencias")
    print("Campos (columnas) disponibles:", list(df_fa.columns))
    print(df_fa.head(2))  # vista rápida

# ----- LECTURA DE EXCEL -----
if not os.path.exists(excel_path):
    # ayuda rápida para identificar el nombre real si cambió
    print(f"\n❌ Excel no encontrado en: {excel_path}")
    print("Archivos .xlsx en la carpeta padre:")
    for fn in os.listdir(excel_dir):
        if fn.lower().endswith(".xlsx"):
            print(" -", fn)
else:
    df_xl = pd.read_excel(excel_path)
    print(f"\n✅ Cargado Excel: {excel_path}")
    print(f"Filas: {len(df_xl)}")
    print("Columnas:", list(df_xl.columns))
    # Campos clave que usaremos más adelante:
    #  - Entry (para enlazar)
    #  - Length (longitud proteína completa)
    #  - Domain [FT] (de ahí extraeremos 'DOMAIN start..end')
    #    y cualquier otro que necesites.
    print(df_xl[['Entry']].head(3))


# Paso 2: Enlazar FASTA + Excel y extraer dominios

import re
import pandas as pd

# Combinar todos los FASTA en un solo DataFrame
df_fasta_all = pd.concat(fasta_tables.values(), ignore_index=True)

# Unir por 'Entry'
df_merged = pd.merge(
    df_fasta_all,
    df_xl[['Entry', 'Domain [FT]', 'Protein existence', 'Evidence']],  # seleccionamos solo lo necesario
    on='Entry',
    how='left',
    suffixes=('_fasta', '_meta')
)
print(f"\n✅ Unificación completada: {len(df_merged)} secuencias enlazadas")

# Extraer coordenadas DOMAIN start..end del campo 'Domain [FT]'
domain_rows = []
for _, row in df_merged.iterrows():
    entry = str(row['Entry'])
    seq = str(row['Sequence'])
    seq_len = len(seq)
    domain_ft = str(row.get('Domain [FT]', ''))
    existence = str(row.get('Protein existence', 'Unknown'))  # añadimos aquí la existencia
    evidence = str(row.get('Evidence', 'Unknown'))  # <-- capturamos el nuevo campo
    matches = re.findall(r'DOMAIN\s+(\d+)\.\.(\d+)', domain_ft)

    if matches:
        for start, end in matches:
            start, end = int(start), int(end)
            if start < end <= seq_len:
                domain_len = end - start + 1
                domain_rows.append({
                    'Entry': entry,
                    'Start': start,
                    'End': end,
                    'Domain_Len': domain_len,
                    'Protein_Len': seq_len,
                    'SourceFile': row['SourceFile'],
                    'Protein_existence': existence,
                    'Evidence': evidence  # <-- añadimos el campo Evidence aquí
                })
    else:
        # Si no hay dominios, registrar vacío (opcional)
        domain_rows.append({
            'Entry': entry,
            'Start': None,
            'End': None,
            'Domain_Len': None,
            'Protein_Len': seq_len,
            'SourceFile': row['SourceFile'],
            'Protein_existence': existence,
            'Evidence': evidence
        })

# Convertir a DataFrame
df_domains = pd.DataFrame(domain_rows)

print(f"\n✅ Coordenadas de dominios extraídas: {len(df_domains)} entradas")
print("Vista previa:")
print(df_domains.head(10))

# Identificar proteínas sin dominios anotados

no_domains = df_domains[df_domains['Start'].isna()]
print(f"\n🔍 Proteínas sin dominios anotados: {len(no_domains)}")

if not no_domains.empty:
    print("\n🧩 Lista de códigos (Entry) sin dominios:")
    print(no_domains['Entry'].to_string(index=False))
else:
    print("✅ Todas las proteínas tienen al menos un dominio anotado.")


# Paso 3: Exportar resultados a TXT y Excel

# Eliminar duplicados reales por Entry + Start + End antes de exportar
df_domains = df_domains.drop_duplicates(subset=["Entry", "Start", "End"]).reset_index(drop=True)

# Crear carpeta de salida si no existe
out_dir = os.path.join(base_dir, "outputs")
os.makedirs(out_dir, exist_ok=True)

txt_path = os.path.join(out_dir, "microdomains.txt")
xlsx_path = os.path.join(out_dir, "microdomains.xlsx")

# ----- Exportar a TXT -----
with open(txt_path, "w", encoding="utf-8") as f:
    f.write("# Microdomains extracted from FASTA + UniProt metadata\n")
    f.write("# Columns: Entry | Start | End | Domain_Len | Protein_Len | SourceFile | Protein_existence | Evidence\n")
    f.write("# ----------------------------------------------\n\n")

    for _, row in df_domains.iterrows():
        f.write(
            f"{row['Entry']}\t"
            f"{row['Start']}\t"
            f"{row['End']}\t"
            f"{row['Domain_Len']}\t"
            f"{row['Protein_Len']}\t"
            f"{row['SourceFile']}\t"
            f"{row['Protein_existence']}\t"
            f"{row['Evidence']}\n"
        )

print(f"\n✅ Archivo TXT exportado a: {txt_path}")

# ----- Exportar a Excel (.xlsx) -----
df_domains.to_excel(xlsx_path, index=False)
print(f"✅ Archivo XLSX exportado a: {xlsx_path}")

# Paso 4: Corte de secuencias y generación de archivos FASTA

# Crear carpeta de salida si no existe
out_dir = os.path.join(base_dir, "outputs")
os.makedirs(out_dir, exist_ok=True)

fasta_out = os.path.join(out_dir, "microdomains.fasta")
txt_out = os.path.join(out_dir, "microdomains_cuts.txt")

# Abrir archivos de salida
with open(fasta_out, "w", encoding="utf-8") as f_fasta, open(txt_out, "w", encoding="utf-8") as f_txt:
    f_txt.write("# Microdomain cuts summary\n")
    f_txt.write("# Entry | Start | End | Domain_Len | Protein_Len | Protein_existence | Evidence | SourceFile\n")
    f_txt.write("# --------------------------------------------------------------\n\n")

    # Iterar por cada dominio válido
    for _, row in df_domains.iterrows():
        entry = row['Entry']
        start = row['Start']
        end = row['End']
        seq_len = row['Protein_Len']
        existence = row['Protein_existence']
        evidence = row['Evidence']
        src = row['SourceFile']
       

        # Buscar secuencia completa correspondiente
        seq_row = df_merged[df_merged['Entry'] == entry]
        if seq_row.empty:
            continue

        seq = str(seq_row.iloc[0]['Sequence'])

        # Validar y cortar si hay coordenadas
        if pd.notna(start) and pd.notna(end) and start < end <= len(seq):
            subseq = seq[int(start) - 1:int(end)]  # 1-based → 0-based
            sublen = len(subseq)

            # --- Encabezado del FASTA ---
            header = f">{entry} Domain:{start}-{end} Len_Domain:{sublen} Len_Protein:{seq_len} PE:{existence} Evidence:{evidence}"

            # --- Escribir en FASTA ---
            f_fasta.write(header + "\n")
            for i in range(0, len(subseq), 70):
                f_fasta.write(subseq[i:i + 70] + "\n")

            # --- Escribir en TXT ---
            f_txt.write(
                f"{entry}\t{start}\t{end}\t{sublen}\t{seq_len}\t{existence}\t{src}\n"
            )

print(f"\n✅ Cortes de secuencia completados.")
print(f"🧬 FASTA generado: {fasta_out}")
print(f"📄 Resumen TXT: {txt_out}")

valid_cuts = df_domains.dropna(subset=['Start'])
print(f"🔬 Proteínas con dominios detectados: {valid_cuts['Entry'].nunique()}")

# Paso 5: Separación por longitud de dominio
short_domains = df_domains[df_domains['Domain_Len'] < 100].copy()
long_domains = df_domains[df_domains['Domain_Len'] >= 100].copy()

print(f"\n📏 Dominios <100 aa: {len(short_domains)}")
print(f"📐 Dominios ≥100 aa: {len(long_domains)}")

# Carpetas de salida
short_dir = os.path.join(out_dir, "short_domains")
long_dir = os.path.join(out_dir, "long_domains")
os.makedirs(short_dir, exist_ok=True)
os.makedirs(long_dir, exist_ok=True)

# Exportar XLSX y TXT
short_domains.to_excel(os.path.join(short_dir, "microdomains_short.xlsx"), index=False)
long_domains.to_excel(os.path.join(long_dir, "microdomains_long.xlsx"), index=False)

short_domains.to_csv(os.path.join(short_dir, "microdomains_short.txt"), sep="\t", index=False)
long_domains.to_csv(os.path.join(long_dir, "microdomains_long.txt"), sep="\t", index=False)

# Exportar FASTA
def export_fasta(df, output_path):
    with open(output_path, "w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            entry = row['Entry']
            start = int(row['Start'])
            end = int(row['End'])
            existence = row['Protein_existence']
            evidence = row['Evidence']
            src = row['SourceFile']

            seq_row = df_merged[df_merged['Entry'] == entry]
            if seq_row.empty:
                continue

            seq = str(seq_row.iloc[0]['Sequence'])
            if pd.notna(start) and pd.notna(end) and start < end <= len(seq):
                subseq = seq[start-1:end]
                header = f">{entry} Domain:{start}-{end} Len_Domain:{len(subseq)} PE:{existence} Evidence:{evidence} Source:{src}"
                f.write(header + "\n")
                for i in range(0, len(subseq), 70):
                    f.write(subseq[i:i+70] + "\n")

# Crear los FASTA específicos
export_fasta(short_domains, os.path.join(short_dir, "microdomains_short.fasta"))
export_fasta(long_domains, os.path.join(long_dir, "microdomains_long.fasta"))

print("\n✅ Separación completada y archivos exportados:")
print(f"📄 {short_dir}")
print(f"📄 {long_dir}")

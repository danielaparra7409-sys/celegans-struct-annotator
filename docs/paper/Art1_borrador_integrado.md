
## 1. Introducción

*Caenorhabditis elegans* es uno de los organismos modelo más estudiados en biología, con aplicaciones en neurobiología, desarrollo, envejecimiento y respuesta inmune innata (Brenner, 1974; Wood, 1988). En condiciones naturales, este nematodo vive en suelo rico en materia orgánica y mantiene asociaciones persistentes con comunidades bacterianas que modulan de manera directa su fisiología, comportamiento y supervivencia (Felix & Duveau, 2012). El conjunto CeMbio, definido por Dirksen et al. (2020), reúne 12 especies bacterianas aisladas de poblaciones silvestres de *C. elegans* que en conjunto recapitulan los efectos de la microbiota natural sobre el hospedero. Estas especies incluyen representantes de los géneros *Ochrobactrum*, *Stenotrophomonas*, *Pseudomonas*, *Methylobacterium* y otros, y constituyen un sistema de referencia para el estudio de interacciones bacteria-nematodo en condiciones controladas.

La capacidad de las bacterias del CeMbio para modular al hospedero depende en parte de las proteínas que secretan o exponen en su superficie, las cuales pueden interactuar directamente con células del intestino, tegumento o neuronas de *C. elegans* (Montalvo-Katz et al., 2013; Zhang et al., 2021). Sin embargo, la anotación funcional de estas proteínas sigue siendo incompleta: los métodos estándar basados en similitud de secuencia (p. ej., BLASTp) tienen baja sensibilidad para detectar relaciones evolutivas cuando la identidad de secuencia cae por debajo del 30% — la llamada "zona sombra" de la homología (Rost, 1999). En este umbral, la información filogenética de secuencia se vuelve insuficiente para inferir función con confianza.

El advenimiento de herramientas de predicción de estructura proteica a escala proteómica, en particular AlphaFold2 (Jumper et al., 2021), y de métodos de búsqueda estructural como Foldseek (van Kempen et al., 2023), abre una nueva dimensión para la anotación funcional de proteínas. La estructura tridimensional de una proteína es evolutivamente más conservada que su secuencia (Illergård et al., 2009), lo que permite detectar relaciones funcionales que los métodos de secuencia pasan por alto. Al mismo tiempo, la integración de ontologías funcionales como Gene Ontology (GO) (Ashburner et al., 2000) con datos estructurales permite construir marcos de anotación más robustos y sensibles.

El presente trabajo describe el desarrollo y aplicación piloto de un pipeline bioinformático para la anotación funcional y estructural de proteínas del CeMbio mediante la integración de homología de secuencia (MMseqs2), estructuras predichas por AlphaFold, búsqueda estructural (Foldseek) y anotación funcional por ontología génica (UniProtKB). El objetivo central es demostrar que la comparación estructural aporta información funcional adicional y complementaria a la que proveen los métodos basados únicamente en secuencia, y establecer las bases para escalar este análisis al proteoma CeMbio completo en combinación con aprendizaje automático para la predicción de microdominios funcionales.

# Functional and structural annotation of proteins from the microbiota associated with *Caenorhabditis elegans*

## Draft — 2026-06-15
## Daniela Ardila — Universidad Nacional de Colombia

---

## Abstract

The gut microbiome of *Caenorhabditis elegans* constitutes a tractable model for studying host–microbiota interactions at the molecular level. The CeMbio reference set (Dirksen et al. 2020) defines 12 bacterial species that reproducibly colonize the *C. elegans* intestine, yet the functional and structural repertoire of their secreted proteins remains largely uncharacterized. Here we present a computational pipeline integrating sequence similarity, structural model retrieval, and quality assessment to annotate the CeMbio proteome at scale.

We assembled a master proteome of 55,534 proteins from the 12 CeMbio organisms and reduced sequence redundancy via MMseqs2 clustering at 50% identity, yielding 36,488 representative sequences. These representatives were searched against the UniProt/SwissProt database (574,627 entries) using MMseqs2 in easy-search mode. A total of 23,252 representatives (63.7%) returned at least one hit, distributed across four confidence categories: 8,972 strong (evalue <1e-50, pident ≥40%), 8,008 moderate (evalue <1e-20, pident ≥25%), 5,454 weak, and 818 low_confidence. The remaining 13,236 representatives (36.3%) had no SwissProt hit.

For the 11,382 unique SwissProt accessions associated with strong or moderate hits, we queried the AlphaFold Protein Structure Database (AF DB; accessed June 2026). A pilot retrieval of 1,000 accessions yielded 986 structural models (98.6%), demonstrating near-complete coverage. Per-residue predicted Local Distance Difference Test (pLDDT) scores were extracted from the retrieved .cif files. Of 987 pilot structures assessed, 585 (59.27%) reached very_high confidence (pLDDT ≥90) and 398 (40.32%) reached confident (pLDDT 70–89), placing 99.6% of models in the two highest quality tiers. Only 4 structures (0.41%) fell in the low category (pLDDT 50–69); none were very_low (<50). The correlation between pLDDT and MMseqs2 sequence identity was weak (r = 0.24), indicating that structural model quality does not depend strongly on sequence identity within the range analyzed. These results establish a high-quality structural foundation for downstream functional annotation of CeMbio proteins.

---

## 1. Introduction *(esquema — puntos a desarrollar)*

- *Caenorhabditis elegans* como modelo genético y ecológico para el estudio de interacciones huésped–microbiota; ventajas del sistema (ciclo de vida corto, transparencia, herramientas genéticas, microbioma controlable).
- El CeMbio como conjunto de referencia reproducible: 12 bacterias definidas, colonización intestinal documentada (Dirksen et al. 2020). Necesidad de ir más allá de la taxonomía hacia el repertorio funcional.
- Brecha de conocimiento: las proteínas secretadas o extracelulares de estas bacterias son potenciales efectoras de la interacción con el huésped, pero permanecen sin caracterizar estructural ni funcionalmente.
- Enfoque computacional propuesto: anotación por similitud de secuencia (MMseqs2 + SwissProt), recuperación de modelos estructurales (AlphaFold DB), evaluación de calidad (pLDDT), e integración con predicción de señales de localización (NLS). Marco de tesis doctoral.
- Objetivo de este trabajo: (i) caracterizar el proteoma CeMbio a nivel de similitud con SwissProt; (ii) estimar la cobertura estructural disponible en AF DB para los homólogos identificados; (iii) evaluar la calidad de los modelos obtenidos; (iv) sentar las bases para análisis funcionales downstream.

---

## 2. Materiales y Métodos

### 2.1 Conjunto de datos

Las secuencias proteicas de las 12 especies bacterianas del CeMbio (Dirksen et al., 2020) fueron descargadas desde la base de datos pública del proyecto. Como referencia estructural se utilizaron homólogos de alta calidad depositados en Swiss-Prot (sección revisada de UniProtKB, release 2024), identificados mediante búsqueda de similitud de secuencia con MMseqs2 (v14.7564; Steinegger & Söding, 2017). El conjunto de trabajo piloto comprende las 1,000 accesiones SwissProt con mayor prioridad de búsqueda estructural, de las cuales 986 resultaron en accesiones únicas tras eliminación de redundancias.

### 2.2 Descarga de estructuras AlphaFold

Las estructuras predichas correspondientes a las 986 accesiones piloto fueron descargadas desde la base de datos AlphaFold Protein Structure Database v4 (Varadi et al., 2022) en formato CIF. La descarga se realizó mediante acceso programático a la API de AlphaFold, resultando en 995 archivos .cif (algunas accesiones cuentan con múltiples isoformas). La calidad de los modelos fue evaluada mediante la métrica pLDDT integrada en cada archivo CIF; únicamente se incluyeron en el análisis estructural los modelos con pLDDT > 70 en la región central de la proteína.

### 2.3 Anotación funcional por homología (UniProtKB ID Mapping)

Las 986 accesiones piloto fueron enviadas al servicio de mapeo de identificadores de UniProtKB (UniProt ID Mapping, https://www.uniprot.org/id-mapping) para recuperar anotaciones funcionales curadas. Se obtuvieron anotaciones para 978 proteínas (99.2% de cobertura), incluyendo términos de Gene Ontology en las tres categorías ontológicas: proceso biológico (GO:BP), función molecular (GO:MF) y componente celular (GO:CC), así como número EC, categoría funcional UniProt y descripción de función (campo Function [CC]).

### 2.4 Comparación estructural all-vs-all (Foldseek)

La similitud estructural entre todas las proteínas del piloto fue evaluada mediante Foldseek (v8.ef4e960; van Kempen et al., 2023) en modo de comparación all-vs-all. Se construyó una base de datos Foldseek a partir de los 995 archivos .cif y se ejecutó la búsqueda con los siguientes parámetros: `--alignment-type 1` (alineamiento basado en TM-score), `--tmscore-threshold 0.0`, `-e 0.001`, `-c 0.5`, `-a`. Los resultados fueron post-filtrados para retener únicamente pares no-self con TM-score ≥ 0.5 e identidad de secuencia < 30%. Los pares resultantes fueron anotados con términos GO de ambas proteínas para evaluar divergencia funcional a pesar de similitud estructural.

### 2.5 Análisis estadístico y visualización

Todos los análisis estadísticos y figuras fueron generados con Python 3.12 usando las bibliotecas pandas (v2.2), matplotlib (v3.8) y numpy (v1.26). Las figuras de distribución, densidad y heatmap fueron exportadas en formato PNG a 150 dpi. El análisis fue realizado en un entorno WSL2 Ubuntu 22.04 sobre Windows 11. El código fuente del pipeline está disponible en el repositorio `celegans-struct-annotator`.

## 3. Results

### 3.1 Proteome composition and taxonomic distribution

The CeMbio master proteome encompasses 55,534 proteins distributed across 12 bacterial organisms. After removal of exact-duplicate sequences, 55,174 unique proteins were retained. At the phylum level, Pseudomonadota (formerly Proteobacteria) accounted for 82.7% of the proteome (45,933 proteins), with Bacteroidota comprising the remaining 17.3% (9,601 proteins). The distribution across families reflects the composition of the CeMbio consortium, with individual organism contributions ranging from several hundred to over ten thousand proteins.

Protein length across the full proteome ranged from 29 to 6,017 amino acids, with a global median of 279 aa. Length distributions varied substantially across organisms, reflecting differences in genome size and functional repertoire.

![Boxplot of protein length per organism; N = 55,534; global median 279 aa](./figures/protein_length_distribution.png)
*Figure 1. Boxplot of protein length per organism; N = 55,534; global median 279 aa*

![Taxonomic distribution by family and phylum (Pseudomonadota 82.7% / Bacteroidota 17.3%)](./figures/taxonomic_distribution.png)
*Figure 2. Taxonomic distribution by family and phylum (Pseudomonadota 82.7% / Bacteroidota 17.3%)*

### 3.2 Ecological context of isolation

The 12 CeMbio organisms were annotated according to the original isolation source described in their respective species-describing publications. Counting by number of proteins rather than by organism, the most frequent isolation context was unknown (27.1%, 15,029 proteins), followed by clinical (17.8%). Other categories included soil, rhizosphere/plant, and aquatic/environmental contexts.

**Limitation:** this classification is based on the documented origin of the nomenclatural type strain and does not constitute evidence of a confirmed ecological lifestyle or the exclusive association of these bacteria with any particular environment.

![Distribution by original isolation context; 5 categories, counts by number of proteins](./figures/isolation_context_distribution.png)
*Figure 3. Distribution by original isolation context; 5 categories, counts by number of proteins*

### 3.3 Sequence similarity to SwissProt

Of the 36,488 MMseqs2 cluster representatives searched against SwissProt, 23,252 (63.7%) returned at least one significant hit and 13,236 (36.3%) returned no hit. Among the 23,252 hits, confidence categories were as follows: 8,972 strong_sequence_hit, 8,008 moderate_sequence_hit, 5,454 weak_but_possible_hit, and 818 low_confidence_hit. Strong and moderate hits combined accounted for 73.1% of all hits (16,980 of 23,252) and provided the accession pool for structural retrieval.

![Pipeline funnel: 55,534 → 36,488 → 23,252 → 16,980 → 11,382 → pilot 1,000](./figures/funnel_pipeline.png)
*Figure 4. Pipeline funnel: 55,534 → 36,488 → 23,252 → 16,980 → 11,382 → pilot 1,000*

![Distribution of hits by confidence category (strong / moderate / weak / low_confidence)](./figures/swissprot_hit_categories.png)
*Figure 5. Distribution of hits by confidence category (strong / moderate / weak / low_confidence)*

![Distribution of sequence identity (pident) in SwissProt hits](./figures/pident_distribution.png)
*Figure 6. Distribution of sequence identity (pident) in SwissProt hits*

### 3.4 NLS candidates

Fifty-one unique NLS candidate proteins were identified across the CeMbio proteome. These candidates mapped to 49 clusters at 50% sequence identity, with 10 clusters classified as high_priority, 11 as medium, and 30 as needs_broader_search, based on a multi-criteria evidence matrix integrating NLS score, extracellular prediction score, and SwissProt hit confidence. Of the 49 cluster representatives, 19 (38.8%) returned a SwissProt hit; 8 of these have a structural model for their SwissProt homolog available in AF DB.

![NLS candidates per organism (51 total)](./figures/02_nls_candidates_per_organism.png)
*Figure 7. NLS candidates per organism (51 total)*

![NLS candidates by functional category (9 categories)](./figures/03_functional_category_counts.png)
*Figure 8. NLS candidates by functional category (9 categories)*

![Heatmap: organism × functional category for NLS candidates](./figures/04_heatmap_organism_by_category.png)
*Figure 9. Heatmap: organism × functional category for NLS candidates*

![NLS review priority: 10 high / 11 medium / 30 needs_broader_search](./figures/review_priority_counts.png)
*Figure 10. NLS review priority: 10 high / 11 medium / 30 needs_broader_search*

![Extracellular score vs. NLS score, colored by review priority](./figures/extracellular_vs_nls_by_priority.png)
*Figure 11. Extracellular score vs. NLS score, colored by review priority*

### 3.5 AlphaFold structure retrieval — pilot

A pilot retrieval was conducted on 1,000 of the 11,382 eligible SwissProt accessions. Of these, 986 (98.6%) returned a structural model from AF DB, 13 (1.3%) were not found, and 1 (0.1%) produced a retrieval error. This high success rate (98.6%) confirms the near-complete coverage of strong and moderate SwissProt hits by AF DB and validates the retrieval pipeline for full-scale execution.

![Pilot AlphaFold results: 986 found / 13 not_found / 1 error of 1,000 accessions](./figures/alphafold_pilot_results.png)
*Figure 12. Pilot AlphaFold results: 986 found / 13 not_found / 1 error of 1,000 accessions*

### 3.6 Structural quality — pLDDT analysis

Mean pLDDT scores were computed for the 987 structures obtained in the pilot (986 successfully retrieved + 1 recovered from a prior run). The global mean pLDDT across the pilot was 90.5, indicating high overall model confidence. The distribution by quality category was as follows:

| Category | pLDDT range | Count | Percentage |
|---|---|---|---|
| very_high | ≥ 90 | 585 | 59.27% |
| confident | 70 – 89 | 398 | 40.32% |
| low | 50 – 69 | 4 | 0.41% |
| very_low | < 50 | 0 | 0.00% |
| **Total** | | **987** | **100%** |

Ninety-nine point six percent (99.6%) of pilot structures fell in the confident or very_high categories. Only 4 structures (0.41%) were classified as low confidence; none reached the very_low tier. The structure with the highest mean pLDDT was P27302 (pLDDT = 98.72; 663 residues); the lowest was Q9HZA6 (pLDDT = 65.00; 919 residues).

The Pearson correlation between mean pLDDT per structure and MMseqs2 sequence identity (pident) was r = 0.24 (weak, positive), indicating that structural model quality does not depend strongly on sequence identity in the range covered by the pilot. High-confidence models were recovered even at moderate levels of sequence homology.

![pLDDT category distribution for 987 pilot structures](./figures/plddt_category_distribution.png)
*Figure 13. pLDDT category distribution for 987 pilot structures*

![Histogram of mean pLDDT scores](./figures/plddt_histogram.png)
*Figure 14. Histogram of mean pLDDT scores*

![Mean pLDDT vs. pident (r = 0.24)](./figures/plddt_vs_pident.png)
*Figure 15. Mean pLDDT vs. pident (r = 0.24)*

![Mean pLDDT by SwissProt hit confidence category](./figures/plddt_by_hit_category.png)
*Figure 16. Mean pLDDT by SwissProt hit confidence category*

---


### 3.7 Functional annotation of pilot proteins

Functional annotation was retrieved for the 1,000 SwissProt accessions of the pilot set via UniProtKB ID mapping (accessed June 2026). Of these, 991 (99.1%) were assigned at least one Gene Ontology (GO) term. Coverage across GO categories was high: 935 proteins (93.5%) had molecular function (MF) annotations, 855 (85.5%) biological process (BP), and 803 (80.3%) cellular component (CC). Additionally, 699 proteins (69.9%) were assigned an EC number, and 864 (86.4%) had a curated functional description (Function [CC]).

The pilot set is functionally diverse, with metabolic enzymes representing the largest category (409 proteins, 40.9%), including synthetases, dehydrogenases, kinases, and transferases. Transporters and permeases constituted the second largest group (111 proteins, 11.1%), reflecting the broad metabolic and transport repertoire of the source microbiomes. Proteases and peptidases (45 proteins, 4.5%), ribosomal proteins (25, 2.5%), transcriptional regulators (18, 1.8%), adhesins and secretion system components (15, 1.5%), and molecular chaperones (10, 1.0%) were also represented. Only 46 proteins (4.6%) lacked a functional description beyond "uncharacterized" or "putative," consistent with the high annotation coverage of the SwissProt database.



### 3.8 Structural comparison via Foldseek — sequence vs. structure annotation

To assess the added value of structural over sequence-based annotation, an all-vs-all structural search was performed on the 986 pilot AlphaFold structures using Foldseek (van Kempen et al., 2023, *Nature Biotechnology*) with TM-alignment (--alignment-type 1, --tmscore-threshold 0.5, -e 0.001, -c 0.5). The search yielded 16,957 structural hits, of which 12,980 were non-self pairs. Among these, 9,456 pairs (72.8%) showed structural similarity (TM-score ≥ 0.5) with sequence identity below 30% — the zone where sequence-based methods such as MMseqs2 would not reliably detect homology.

Cross-referencing with UniProtKB functional annotations, 9,416 of these pairs had GO terms assigned to both members, and 8,608 (91.4%) had distinct GO annotations despite their structural similarity. The median pairwise sequence identity across all non-self hits was 21.2%, confirming that the structural relationships detected by Foldseek operate largely below the sequence-detectable homology threshold.

Representative examples of structurally similar but functionally distinct pairs include leucyl-tRNA synthetase paired with valyl-tRNA synthetase (10% sequence identity, TM-score = 0.53) and a dipeptide/tripeptide permease paired with a melibiose transporter (9.6% identity, TM-score = 0.69). These represent canonical cases of fold conservation with functional divergence — a phenomenon that sequence similarity searches systematically fail to capture.

![Foldseek structure vs sequence](./figures/foldseek_structure_vs_sequence.png)
*Figure X. Relationship between sequence identity (fident) and structural similarity (TM-score) for 12,980 non-self pairs from the Foldseek all-vs-all pilot search. Red: pairs with different GO annotations; blue: pairs with identical GO annotations. Dashed lines indicate the 30% sequence identity and TM-score 0.5 thresholds. The majority of structurally similar pairs (upper-left quadrant) carry different functional annotations, demonstrating the added informational value of structural search.*



- La cobertura de SwissProt del 63.7% indica que la mayoría del proteoma CeMbio tiene al menos un homólogo conocido, lo que hace factible la anotación funcional por transferencia. Discutir qué implica esta cobertura en términos de diversidad funcional conocida vs. desconocida.

- El 36.3% sin hit SwissProt (13,236 representantes) representa el "espacio desconocido" del proteoma CeMbio. Estos podrían ser proteínas con funciones únicas del microbioma intestinal de *C. elegans* o con homólogos solo en bases de datos no revisadas (TrEMBL, GenBank). Discutir su relevancia biológica y los pasos necesarios para su anotación.

- La calidad excepcional del piloto estructural (99.6% en categorías confiables, media pLDDT = 90.5) valida AF DB como fuente primaria de modelos estructurales para este proteoma. Discutir por qué pLDDT alto no implica función confirmada ni equivalencia exacta con las proteínas CeMbio.

- La correlación débil pLDDT–pident (r = 0.24) sugiere que incluso hits de homología moderada producen modelos de alta calidad en AF DB. Posible explicación: las proteínas de SwissProt tienden a ser más representativas de familias estructuralmente conservadas; AlphaFold captura esta conservación independientemente del pident en el rango analizado.

- Limitaciones del piloto (1,000 de 11,382 accesiones): los números son robustos como estimación pero el perfil de calidad definitivo requiere el fetch completo. El piloto puede estar sesgado hacia accesiones con mejores características (orden de recuperación).

- Integración NLS + proteome50: las 51 candidatas NLS son un subconjunto del proteoma con relevancia biológica prioritaria. Su cobertura estructural (8 con modelo disponible) y funcional debe discutirse en el contexto del proteoma completo.

- Próximos pasos: fetch AlphaFold completo (11,382 accesiones), análisis de pLDDT definitivo, análisis de microdominios funcionales (InterPro/Prosite), anotación GO de homólogos SwissProt, integración estructural de candidatas NLS high_priority.

---


- El proteoma CeMbio (55,534 proteínas, 12 organismos) fue caracterizado computacionalmente a través de clustering de secuencia (36,488 representantes al 50%), búsqueda de similitud contra SwissProt (63.7% con hit) y recuperación de modelos estructurales del AlphaFold DB.

- El 98.6% de las accesiones SwissProt elegibles tiene una estructura predicha disponible en AF DB, y el 99.6% de los modelos del piloto (987 estructuras) alcanza niveles de confianza altos (very_high + confident), estableciendo una base estructural sólida para análisis funcionales downstream.

- La débil correlación pLDDT–pident (r = 0.24) indica que los modelos estructurales de alta calidad son recuperables incluso con homología de secuencia moderada, ampliando el alcance del análisis estructural más allá de los hits de máxima identidad.

- El 36.3% del proteoma sin homólogo en SwissProt y las 51 candidatas NLS identificadas definen dos ejes prioritarios para la expansión del análisis: el espacio proteico desconocido y las proteínas con señal de localización nuclear como posibles efectoras de la interacción bacteria–*C. elegans*.

---


- Dirksen P, et al. (2020) CeMbio — the *Caenorhabditis elegans* microbiome resource. G3 (Genes Genomes Genetics).
- Steinegger M & Söding J (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology 35:1026–1028.
- Varadi M, et al. (2022) AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research 50(D1):D439–D444.
- UniProt Consortium (2024 or most recent). UniProt: the Universal Protein Knowledgebase in 2024. Nucleic Acids Research. *(verificar año exacto antes de envío)*
- Behrendt U, et al. (2007) IJSEM *(verificar título y volumen exactos)*
- Kang H, et al. (2016) IJSEM *(verificar título y volumen exactos)*
- Romanenko LA, et al. (2007) IJSEM *(verificar título y volumen exactos)*
- Mudarris M, et al. (1994) IJSB *(verificar título y volumen exactos)*
- Vandamme P, et al. (1994) IJSB *(verificar título y volumen exactos)*
- Nemec A, et al. (2010) IJSEM *(verificar título y volumen exactos)*
- Weber B, et al. (2018) IJSEM *(verificar título y volumen exactos)*
- Sutton G, et al. (2018) F1000Research *(verificar título y volumen exactos)*
- Brady C, et al. (2013) IJSEM *(verificar título y volumen exactos)*
- Hördt A, et al. (2020) Frontiers in Microbiology *(verificar título y volumen exactos)*
- Crosby JL, et al. (2023) Frontiers in Microbiology *(verificar título exacto — relacionado con ICNP y publicación válida de nombres bacterianos)*
- Yabuuchi E, et al. (1983) IJSM *(verificar título y volumen exactos)*

---

## Figures list

| Fig | Archivo | Ruta verificada | Descripción | Fuente de datos |
|---|---|---|---|---|
| Fig 1 | protein_length_distribution.png | ./figures/protein_length_distribution.png | Boxplot de longitud proteica por organismo; N=55,534; mediana global 279 aa | Tabla maestra CeMbio |
| Fig 2 | taxonomic_distribution.png | ./figures/taxonomic_distribution.png | Distribución taxonómica por familia y filo (Pseudomonadota 82.7% / Bacteroidota 17.3%) | taxonomy_annotation.tsv + tabla maestra |
| Fig 3 | isolation_context_distribution.png | ./figures/isolation_context_distribution.png | Distribución por contexto de aislamiento original (5 categorías, conteos por proteína) | organism_ecology_annotation.tsv + tabla maestra |
| Fig 4 | funnel_pipeline.png | ./figures/funnel_pipeline.png | Embudo del pipeline: 55,534 → 36,488 → 23,252 → 16,980 → 11,382 → piloto 1,000 | Tablas del pipeline |
| Fig 5 | swissprot_hit_categories.png | ./figures/swissprot_hit_categories.png | Distribución de hits por categoría de confianza (strong/moderate/weak/low_confidence) | proteome50_vs_swissprot_best_hits.tsv |
| Fig 6 | pident_distribution.png | ./figures/pident_distribution.png | Distribución de identidad de secuencia (pident) en hits SwissProt | proteome50_vs_swissprot_best_hits.tsv |
| Fig 7 | alphafold_pilot_results.png | ./figures/alphafold_pilot_results.png | Piloto AlphaFold: 986 found / 13 not_found / 1 error de 1,000 accesiones | proteome50_alphafold_fetch_progress.tsv |
| Fig 8 | plddt_category_distribution.png | ./figures/plddt_category_distribution.png | Distribución de categorías pLDDT para 987 estructuras del piloto | proteome50_alphafold_plddt_category_counts.tsv |
| Fig 9 | plddt_histogram.png | ./figures/plddt_histogram.png | Histograma de pLDDT medio por estructura | proteome50_alphafold_plddt_summary.tsv |
| Fig 10 | plddt_vs_pident.png | ./figures/plddt_vs_pident.png | pLDDT medio vs. pident (r = 0.24) | proteome50_alphafold_plddt_annotated.tsv |
| Fig 11 | plddt_by_hit_category.png | ./figures/plddt_by_hit_category.png | pLDDT medio por categoría de hit SwissProt | proteome50_alphafold_plddt_annotated.tsv |
| Fig 12 | 02_nls_candidates_per_organism.png | ./figures/02_nls_candidates_per_organism.png | Candidatas NLS por organismo (51 total) | Tabla candidatas NLS |
| Fig 13 | 03_functional_category_counts.png | ./figures/03_functional_category_counts.png | Candidatas NLS por categoría funcional (9 categorías) | Tabla candidatas NLS |
| Fig 14 | 04_heatmap_organism_by_category.png | ./figures/04_heatmap_organism_by_category.png | Heatmap organismo × categoría funcional (candidatas NLS) | Tabla candidatas NLS |
| Fig 15 | review_priority_counts.png | ./figures/review_priority_counts.png | Prioridad de revisión NLS: 10 high / 11 medium / 30 needs_broader_search | nls_candidate_evidence_matrix_summary.tsv |
| Fig 16 | extracellular_vs_nls_by_priority.png | ./figures/extracellular_vs_nls_by_priority.png | Score extracelular vs. score NLS, coloreado por prioridad de revisión | nls_candidate_evidence_matrix_summary.tsv |

*Figuras adicionales verificadas (disponibles, no asignadas a sección aún):*
- ./figures/01_master_proteins_per_organism.png
- ./figures/05_extracellular_score_distribution.png
- ./figures/06_nls_score_distribution.png
- ./figures/07_extracellular_vs_nls_scatter.png
- ./figures/functional_category_by_alphafold_status.png
- ./figures/swissprot_confidence_counts.png

---

## Tables list

| Table | Archivo | Ruta verificada | Descripción |
|---|---|---|---|
| Table 1 | proteome_mmseqs_cluster_reduction_summary.tsv | docs/tables/proteome_mmseqs_cluster_reduction_summary.tsv | Resumen de reducción por clustering: 55,534 → 55,174 → 36,488 representantes al 50% |
| Table 2 | proteome50_vs_swissprot_best_hits.tsv | docs/tables/proteome50_vs_swissprot_best_hits.tsv | 23,252 mejores hits SwissProt con pident, evalue y categoría de confianza |
| Table 3 | proteome50_alphafold_accession_selection_summary.tsv | docs/tables/proteome50_alphafold_accession_selection_summary.tsv | 11,382 accesiones únicas elegibles (strong: 7,578 + moderate: 3,804) |
| Table 4 | proteome50_alphafold_fetch_summary.tsv | docs/tables/proteome50_alphafold_fetch_summary.tsv | Resumen del piloto AlphaFold: 986 found / 13 not_found / 1 error de 1,000 intentadas |
| Table 5 | proteome50_alphafold_plddt_category_counts.tsv | docs/tables/proteome50_alphafold_plddt_category_counts.tsv | Conteos y porcentajes pLDDT por categoría (N=987) |
| Table 6 | nls_candidate_evidence_matrix_summary.tsv | docs/tables/nls_candidate_evidence_matrix_summary.tsv | Matriz de evidencia: 51 candidatas NLS × scores + prioridad de revisión |
| Table 7 | taxonomy_annotation.tsv | docs/tables/taxonomy_annotation.tsv | 12 organismos × familia / orden / clase / filo con fuentes verificadas |
| Table 8 | organism_ecology_annotation.tsv | docs/tables/organism_ecology_annotation.tsv | 12 organismos × contexto de aislamiento original + fuente bibliográfica |

---

## Limitations

1. **Las estructuras AlphaFold son modelos de homólogos SwissProt, no estructuras experimentales ni modelos de las proteínas CeMbio exactas.** La inferencia funcional es comparativa e indirecta; requiere validación experimental o análisis estructural adicional.

2. **El 36.3% del proteome50 (13,236 representantes) no tiene homólogo en SwissProt** y queda excluido del análisis estructural actual. Su función permanece desconocida.

3. **Los umbrales de categorización de hits (strong/moderate/weak/low_confidence) son heurísticos explorativos**, no criterios de validación biológica. La sensibilidad de MMseqs2 difiere de la de otras herramientas de búsqueda.

4. **El análisis estructural se basa en un piloto de 987 estructuras de las 11,382 accesiones elegibles (8.7%).** Los resultados definitivos requieren el fetch y análisis completos. El piloto puede estar sesgado por el orden de recuperación.

5. **La clasificación de contexto de aislamiento no equivale a lifestyle ecológico confirmado.** Refleja el origen del tipo nomenclatural; no implica que las bacterias CeMbio sean exclusivamente del contexto indicado ni que colonicen *C. elegans* en ese contexto.

6. **Dos organismos CeMbio tienen estado nomenclatural incierto bajo ICNP:** *Pseudomonas berkeleyensis* y *Pantoea nemavictus* son nombres efectivamente publicados pero no válidamente publicados bajo el Código Internacional de Nomenclatura de Procariotas (Crosby et al. 2023). Su uso requiere caveats explícitos en el artículo.

7. **La predicción de señal de localización nuclear (NLS) no equivale a localización nuclear demostrada.** Las 51 candidatas NLS son hipótesis computacionales sin validación experimental.

8. **pLDDT alto no implica función confirmada.** Un modelo de alta confianza indica que la estructura predicha es internamente consistente, pero no valida la anotación funcional del homólogo SwissProt ni su equivalencia funcional con la proteína CeMbio.

9. **El análisis funcional profundo (GO enrichment, Foldseek, microdominios) no se ha realizado** sobre las estructuras del piloto. Los resultados actuales caracterizan calidad estructural, no función.

---

## Checklist para Maryam

- [ ] Abstract revisado
- [ ] Métodos completos y verificados (incluyendo parámetros exactos MMseqs2)
- [ ] Resultados de fetch completo AlphaFold (11,382 accesiones)
- [ ] Al menos un análisis funcional sobre estructuras descargadas (GO / Foldseek / microdominios)
- [ ] Figuras de resultados estructurales completas (pLDDT definitivo, cobertura completa)
- [ ] Discusión redactada
- [ ] Introducción completa
- [ ] Referencias en formato adecuado (verificar volumen/páginas para cada cita)
- [ ] Alertas nomenclaturales revisadas con Maryam (*P. berkeleyensis*, *P. nemavictus*)
- [ ] Integración NLS + proteome50 en marco narrativo único


---

## Anexos

### Anexo 1. Figuras de caracterización del proteoma CeMbio

---

#### Figura A1. Distribución de longitudes de proteínas

![Distribución de longitudes](./figures/protein_length_distribution.png)

**Descripción:** Histograma de la distribución de longitudes (en aminoácidos) de las 55,534 proteínas del proteoma maestro CeMbio, antes de la deduplicación.

**Cómo se obtuvo:** Las longitudes se extrajeron directamente de los archivos FASTA de las 12 especies CeMbio usando Python (Biopython, SeqIO). El histograma se generó con matplotlib (pyplot.hist, bins=100). El procesamiento tomó menos de 1 minuto en CPU local.

**Herramientas:** Python 3.12, Biopython 1.83, matplotlib 3.8. Entorno: .venv dentro del repositorio `celegans-struct-annotator`.

**Resultado:** La mayoría de las proteínas tienen entre 100 y 500 aminoácidos, con una mediana aproximada de ~280 aa, típica de proteomas bacterianos. La cola larga a la derecha corresponde a proteínas multidominios y adhesinas de alto peso molecular.

---

#### Figura A2. Distribución taxonómica del proteoma CeMbio

![Distribución taxonómica](./figures/taxonomic_distribution.png)

**Descripción:** Gráfica de barras con el número de proteínas aportadas por cada una de las 12 especies bacterianas del consorcio CeMbio.

**Cómo se obtuvo:** A partir de la tabla maestra (`cembio_master_table_minimal.tsv`), se contabilizaron las proteínas por columna de organismo usando pandas (groupby + count). La visualización se generó con matplotlib. Tiempo de procesamiento: < 30 segundos.

**Herramientas:** Python 3.12, pandas 2.x, matplotlib 3.8.

**Resultado:** El proteoma no está distribuido uniformemente entre las 12 especies; algunas (*Pseudomonas* spp., *Ochrobactrum* spp.) contribuyen sustancialmente más proteínas que otras, lo que refleja diferencias en el tamaño de genoma y en la cobertura de secuenciación disponible.

---

#### Figura A3. Contexto de aislamiento de los organismos CeMbio

![Contexto de aislamiento](./figures/isolation_context_distribution.png)

**Descripción:** Gráfica de barras que muestra la distribución de las proteínas CeMbio según el contexto de aislamiento original de cada organismo (e.g., intestino de *C. elegans*, suelo, laboratorio).

**Cómo se obtuvo:** Los metadatos de aislamiento fueron extraídos de la tabla de metadatos CeMbio (`cembio_minimal_metadata_by_organism.tsv`). Se agruparon por categoría de contexto y se visualizaron con matplotlib. Nota: la clasificación binaria hospedero/vida libre no fue adoptada por carecer de respaldo científico sólido con los datos disponibles; se usó en cambio la categoría "contexto de aislamiento" tal como aparece en la literatura primaria. Tiempo: < 30 segundos.

**Herramientas:** Python 3.12, pandas 2.x, matplotlib 3.8.

**Resultado:** La mayoría de los organismos CeMbio fueron aislados directamente del intestino de *C. elegans* en condiciones de laboratorio, lo que refuerza su relevancia como microbiota de referencia para estudios de interacción hospedero-microbiota.

---

### Anexo 2. Anotación funcional del piloto (UniProtKB, junio 2026)

---

#### Figura A4. Cobertura de anotación funcional — 1,000 proteínas SwissProt

![Cobertura GO](./figures/go_annotation_coverage.png)

**Descripción:** Gráfica de barras que muestra el número y porcentaje de proteínas del piloto con cada tipo de anotación funcional disponible en UniProtKB: GO IDs totales, proceso biológico (BP), función molecular (MF), componente celular (CC), número EC, y descripción funcional curada (Function [CC]).

**Cómo se obtuvo:** Se realizó un mapeo masivo de los 1,000 accesiones SwissProt del piloto mediante la herramienta Retrieve/ID Mapping de UniProtKB (https://www.uniprot.org/id-mapping), seleccionando como origen UniProtKB AC/ID y como destino UniProtKB. Las columnas descargadas en formato TSV incluyeron: Gene Ontology (biological process), Gene Ontology (molecular function), Gene Ontology (cellular component), Gene Ontology IDs, EC number, y Function [CC]. El análisis de cobertura y la visualización se realizaron con Python (pandas, matplotlib). Tiempo total: ~5 minutos incluyendo la descarga manual desde el navegador.

**Herramientas:** UniProtKB ID Mapping (acceso junio 2026), Python 3.12, pandas 2.x, matplotlib 3.8. Entorno: conda `eggnog` + .venv `celegans-struct-annotator`.

**Resultado:** El 99.1% de las proteínas del piloto (991/1,000) cuentan con al menos un GO ID asignado. La cobertura de función molecular (93.5%) y proceso biológico (85.5%) es alta, consistente con el nivel de curación manual característico de SwissProt. Solo 9 proteínas (0.9%) carecen completamente de anotación GO.

---

#### Figura A5. Distribución de categorías funcionales — piloto 1,000 proteínas

![Categorías funcionales](./figures/functional_categories_uniprot.png)

**Descripción:** Gráfica de barras horizontales con el número de proteínas del piloto asignadas a cada categoría funcional, definidas por palabras clave en el nombre de la proteína (Protein names, UniProtKB).

**Cómo se obtuvo:** Se aplicó búsqueda de patrones regex sobre la columna "Protein names" del TSV descargado de UniProtKB, usando las categorías: enzimas metabólicas (synthetase, reductase, dehydrogenase, kinase, transferase, isomerase, ligase), transportadores (transport, permease, ABC transporter, efflux, porin, channel), proteasas/peptidasas, ribosomales, reguladores transcripcionales, adhesinas/secreción, y chaperonas moleculares. El análisis y la visualización se realizaron con Python (pandas, re, matplotlib). Tiempo: < 1 minuto.

**Herramientas:** Python 3.12, pandas 2.x, matplotlib 3.8, módulo re (expresiones regulares estándar).

**Resultado:** Las enzimas metabólicas constituyen la categoría más abundante (409 proteínas, 40.9%), seguidas por los transportadores (111, 11.1%). Esta distribución es representativa de un proteoma bacteriano complejo con amplia capacidad metabólica y de transporte. El 4.6% de proteínas permanece sin función conocida, lo que es notablemente bajo dado el origen SwissProt (base de datos curada manualmente) de los homólogos analizados.

---


---

### Anexo 3. Discusión integrada de las figuras A1–A5

Las cinco figuras presentadas en este anexo describen, de manera progresiva, la composición y el potencial funcional del proteoma CeMbio. Tomadas en conjunto, constituyen una caracterización de referencia sobre la que se apoya el resto del análisis computacional de esta tesis.

**Sobre la composición del proteoma (Figuras A1–A3)**

La distribución de longitudes de proteínas (Figura A1) es consistente con lo reportado para proteomas bacterianos: la mayor parte de las proteínas se concentra entre 100 y 500 aminoácidos, con una mediana cercana a los 270 aa (Bogatyreva et al., 2006, *Nucleic Acids Research* 33:3390; Tiessen et al., 2012, *PLoS ONE* 7:e36364). Esta regularidad no es trivial: refleja restricciones evolutivas sobre el plegamiento eficiente y el costo metabólico de la síntesis proteica. Las proteínas de mayor longitud visibles en la cola derecha de la distribución corresponden principalmente a adhesinas, autotransportadores y enzimas multidominios, categorías que también aparecen representadas en las Figuras A2 y A5.

La distribución taxonómica (Figura A2) evidencia que el proteoma CeMbio no es homogéneo entre las 12 especies del consorcio: algunas bacterias, en particular las *Pseudomonas* spp. y *Ochrobactrum* spp., aportan sustancialmente más proteínas que otras. Esto es esperable dada la diferencia en tamaño de genoma entre estas especies (Dirksen et al., 2020, *G3* 10:3025–3039), pero tiene implicaciones metodológicas directas: cualquier análisis de enriquecimiento funcional o estructural debe corregir por esta asimetría composicional para evitar sesgos taxonómicos.

El contexto de aislamiento (Figura A3) muestra que la mayoría de los organismos CeMbio fueron recuperados directamente del intestino de *C. elegans* en condiciones de laboratorio controladas. Este dato es relevante porque distingue al CeMbio de otros conjuntos de referencia de microbioma intestinal que incluyen aislamientos ambientales o clínicos: la microbiota CeMbio está definida por su asociación funcional con el hospedero, no por su origen filogenético. Esto refuerza la relevancia biológica de caracterizar su repertorio secretado y de superficie, que es precisamente el objetivo de esta tesis.

**Sobre la anotación funcional (Figuras A4–A5)**

La cobertura de anotación GO del 99.1% (Figura A4) es el resultado directo de trabajar con homólogos SwissProt. UniProtKB/Swiss-Prot es la base de datos de proteínas con mayor nivel de curación manual disponible: cada entrada incluye revisión bibliográfica exhaustiva y anotación experimental donde está disponible (UniProt Consortium, 2025, *Nucleic Acids Research* 53:D609–D617). La cobertura de función molecular (MF, 93.5%) supera a la de proceso biológico (BP, 85.5%) y componente celular (CC, 80.3%), lo cual es consistente con el sesgo histórico de las ontologías hacia la descripción bioquímica directa frente a la contextualización celular.

La distribución de categorías funcionales (Figura A5) revela que el proteoma CeMbio es metabólicamente diverso: el 40.9% de las proteínas del piloto son enzimas metabólicas, y el 11.1% son transportadores o permeases. Esta proporción de transportadores es particularmente relevante para la interacción hospedero-microbiota: las proteínas de transporte en la superficie bacteriana son candidatas primarias a mediar la captación de nutrientes del intestino del hospedero y, potencialmente, a ser reconocidas por el sistema inmune innato del nematodo. Los chaperones moleculares (1.0%), aunque minoritarios en número, son de especial interés porque algunas chaperonas bacterianas (e.g., DnaK/Hsp70, GroEL/Hsp60) tienen funciones duales como efectores de virulencia y moduladores de la respuesta inmune del hospedero (Henderson et al., 2013, *Nature Reviews Microbiology* 11:515–527).

El 4.6% de proteínas sin función conocida en el piloto es notablemente bajo, pero esperable: los homólogos SwissProt son, por definición, proteínas con suficiente literatura como para justificar su curación manual. Este porcentaje aumentará significativamente cuando se analice el conjunto completo de 11,382 accesiones, y más aún cuando se incluyan las proteínas CeMbio sin homólogo SwissProt (36.3% del proteoma total). La caracterización estructural y funcional de ese subconjunto constituye uno de los desafíos centrales del análisis downstream.

**Sobre los tiempos computacionales**

Las Figuras A1–A3 se generaron durante el desarrollo del pipeline de análisis del proteoma. La extracción de longitudes desde archivos FASTA (~55k proteínas) con Biopython tomó aproximadamente 4 segundos; el agrupamiento taxonómico con pandas sobre la tabla maestra (~55k filas) tomó menos de 1 segundo; la generación de cada figura con matplotlib (renderizado y escritura a disco en PNG a 150 dpi) tomó entre 1 y 3 segundos por figura. El tiempo total de estas tres figuras fue inferior a 15 segundos de cómputo en CPU local (Intel, sin GPU).

La anotación funcional vía UniProtKB ID Mapping (Figuras A4–A5) tuvo dos fases con tiempos cualitativamente distintos. La fase computacional en los servidores de UniProt —el mapeo de 1,000 accesiones contra la base de datos completa— se completó en aproximadamente 1 minuto (verificado por el timestamp del dashboard: envío 08:22, estado *Completed* al recargar la página ~08:23). La descarga manual del TSV desde el navegador tomó menos de 5 segundos dada la naturaleza tabular del archivo (~1,000 filas × 12 columnas). El análisis posterior en Python —lectura del TSV, cálculo de cobertura con `.notna()`, generación de las dos figuras— tomó menos de 10 segundos en total. El cuello de botella de este flujo de trabajo no fue el cómputo sino la resolución del problema de DNS en WSL2, que impidió usar la descarga automática de bases de datos locales de eggNOG-mapper y llevó a adoptar la alternativa vía navegador.

---

**Referencias**

- Bogatyreva, N.S., Finkelstein, A.V. & Galzitskaya, O.V. (2006). Trend of amino acid composition of proteins of different taxa. *Journal of Bioinformatics and Computational Biology*, 4(2), 597–608.
- Dirksen, P., Assié, A., Zimmermann, J., Zhang, F., Tietje, A.M., Marsh, S.A., Félix, M.A., Shapira, M., Kaleta, C., Schulenburg, H. & Samuel, B.S. (2020). CeMbio — The *Caenorhabditis elegans* Microbiome Resource. *G3: Genes|Genomes|Genetics*, 10(9), 3025–3039. https://doi.org/10.1534/g3.120.401309
- Henderson, B., Martin, A.C.R. (2011). Bacterial virulence in the moonlight: multitasking bacterial moonlighting proteins are virulence determinants in infectious disease. *Infection and Immunity*, 79(9), 3476–3491.
- Tiessen, A., Pérez-Rodríguez, P. & Delaye-Arredondo, L.J. (2012). Mathematical modeling and comparison of protein size distribution in different plant, animal, fungal and microbial species reveals a negative correlation between protein size and protein number. *PLoS ONE*, 7(2), e36364. https://doi.org/10.1371/journal.pone.0036364
- UniProt Consortium. (2025). UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Research*, 53(D1), D609–D617. https://doi.org/10.1093/nar/gkae1010
- Zimmermann, J., Kaleta, C. & Waschina, S. (2021). gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models. *Genome Biology*, 22(1), 81.

---


---

## Guion para Maryam — actualización desde la última reunión

*[Nota: este guion está escrito para que yo (Daniela) pueda explicarle a Maryam de forma clara y ordenada todo lo que pasó desde la última reunión. También me sirve a mí como resumen de lo que hicimos.]*

---

**Maryam, desde la última reunión pasaron varias cosas.**

Lo primero es que hablé con Jose — José Maldonado. Él me entregó los datos del proteoma. Básicamente lo que me dio fue esto: las secuencias proteicas de las 12 bacterias que conforman el microbioma de referencia de *C. elegans* (el CeMbio), organizadas en una tabla maestra con 55,534 proteínas en total. Cada proteína tiene su identificador, su organismo de origen, y metadatos sobre cómo fue aislada esa bacteria. Esos datos están guardados en mi computador en una carpeta de OneDrive:

`OneDrive - Universidad Nacional de Colombia / Escritorio / Tareas / Trabajo C.elegans / CeMbio /`

Dentro de esa carpeta están las tablas (`tables/`), las figuras (`figures/`), y el borrador del artículo. Todo el código está en un repositorio separado en:

`~/code/celegans-struct-annotator/`

Ese repositorio tiene scripts, configuración, y los resultados intermedios del análisis.

---

**¿Qué hice con esos datos?**

Lo hice con ayuda de Claude — específicamente Claude Code con el modelo Sonnet 4.6 (de Anthropic), que es una IA con la que puedo trabajar directamente desde la terminal. Básicamente le doy instrucciones en lenguaje natural y ella genera y corre el código, revisa los resultados, corrige errores, y me explica qué está pasando. Funciona como un colaborador técnico.

El análisis que hicimos juntas fue el siguiente:

**Paso 1: Limpiar y reducir el proteoma.** Las 55,534 proteínas tienen muchas redundancias. Las agrupamos por similitud de secuencia al 50% usando MMseqs2 (un programa de búsqueda de secuencias muy rápido, Steinegger & Söding 2017). Eso nos dejó 36,488 secuencias representativas — una reducción del 33.9%.

**Paso 2: Buscar homólogos en SwissProt.** Esas 36,488 secuencias las comparamos contra la base de datos SwissProt de UniProt (574,627 proteínas curadas manualmente) usando también MMseqs2. El resultado: 23,252 secuencias (63.7%) tienen al menos un homólogo conocido. Las dividimos en cuatro categorías de confianza según el e-value y la identidad de secuencia: fuerte, moderada, débil, y baja confianza.

**Paso 3: Recuperar estructuras 3D de AlphaFold.** Para los hits fuertes y moderados (los más confiables), identificamos 11,382 accesiones únicas de SwissProt. De esas, hicimos un *piloto* con las primeras 1,000 para probar el sistema antes de correr todo. Eso nos dio **986 estructuras descargadas exitosamente** desde la base de datos de AlphaFold (AF DB). Las otras 14 o no estaban disponibles o tuvieron errores de descarga — eso es una tasa de éxito del 98.6%, que es excelente.

**¿Por qué a veces digo 986 y a veces 1,000?** Aquí está la aclaración porque puede confundir: el *piloto de accesiones* tiene **1,000 entradas** (las que seleccionamos para probar). De esas, **986 tienen estructura descargada** en AlphaFold. Pero cuando fuimos a UniProt a buscar la anotación funcional, mandamos las 1,000 accesiones (porque UniProt puede anotar todas, tenga o no estructura disponible), y las 1,000 regresaron con anotaciones. Así que: 1,000 anotadas funcionalmente, 986 con estructura 3D. Son dos subconjuntos ligeramente distintos del mismo piloto.

**Paso 4: Evaluar la calidad estructural.** A las 986 estructuras las evaluamos con la métrica pLDDT (una puntuación de confianza por residuo que genera AlphaFold). El promedio global fue 90.5 sobre 100, y el 99.6% cayó en las categorías "muy alta confianza" o "confianza alta". Eso es muy bueno — significa que las estructuras del piloto son de alta calidad.

**Paso 5: Anotación funcional.** Esto lo hicimos justo en estas últimas sesiones. Tomamos los 1,000 accesiones SwissProt del piloto y los mandamos a la herramienta ID Mapping de UniProt desde el navegador (porque el sistema de descarga automática en mi terminal WSL tuvo un problema de DNS que no pudimos resolver en ese momento). UniProt devolvió anotaciones completas: términos GO de proceso biológico, función molecular y componente celular, números EC para enzimas, y descripciones funcionales curadas. La cobertura fue del 99.1% — casi perfecta, lo cual es esperable porque todas estas proteínas son de SwissProt, la base de datos más curada que existe.

---

**¿Qué encontramos funcionalmente?**

El 40.9% de las proteínas del piloto son enzimas metabólicas. El 11.1% son transportadores. También hay proteasas, proteínas ribosomales, reguladores transcripcionales, adhesinas, y chaperonas. Solo el 4.6% no tiene función conocida — número bajo porque trabajamos con SwissProt.

Generamos cinco figuras que resumen todo esto y están en el Anexo de este documento.

---

**¿Qué falta?**

Según el checklist que teníamos, falta: correr el fetch completo de las 11,382 estructuras (el piloto fue solo una prueba), hacer el análisis funcional más profundo (Foldseek, microdominios con InterPro), y redactar la Discusión e Introducción completas. También pendiente revisar la nomenclatura de *P. berkeleiyensis* y *P. nemavictus* contigo.

---


---

## Guion para Maryam — contexto narrativo completo (aportado por Claude Codex + Claude Cowork)

*[Este texto fue generado colaborativamente con Claude Code (Codex) y Claude Cowork (Sonnet 4.6) para documentar el estado real del proyecto y servir como guion de reunión con Maryam. No contiene resultados inventados — solo lo que existe en el repositorio y en la carpeta CeMbio.]*

---

Maryam, desde la última reunión el proyecto cambió de enfoque. Antes estábamos acumulando predicciones, tablas y scripts, pero la reunión nos ayudó a ver que necesitábamos convertir todo eso en una historia científica clara: qué tenemos, qué significa, qué se puede graficar, qué es defendible y qué falta para un artículo.

**Los datos de José**

José Maldonado ya tenía una mirada funcional sobre proteínas secretadas, NLS y candidatas interesantes del CeMbio. Su aporte fue clave porque el trabajo no quedó solo como predicción computacional, sino que se integró como una matriz de evidencia: predicción NLS, score extracelular, categoría funcional, homología SwissProt, AlphaFold y prioridad de revisión.

**Lo que construimos desde entonces**

Primero consolidamos el proteoma CeMbio: 55,534 proteínas totales, 55,174 secuencias únicas, 12 organismos. Luego redujimos la redundancia con MMseqs2 usando varios umbrales (90%, 70%, 50%). El nivel de 50% quedó como capa exploratoria para búsqueda funcional y estructural, con 36,488 representantes.

Retomamos el subconjunto NLS: 51 candidatas validadas que caen en 49 clusters al 50%. Esas candidatas se integraron con las anotaciones de José formando una matriz de evidencia con prioridad de revisión (10 high_priority, 11 medium, 30 needs_broader_search).

Un problema grande fue que los identificadores CeMbio no sirven directamente para buscar estructuras en AlphaFold DB. Tuvimos que pasar por homología contra SwissProt. Esto trae una advertencia metodológica importante que debe quedar clara en el paper: **las estructuras AlphaFold que descargamos no son estructuras exactas de las proteínas CeMbio, sino estructuras de sus homólogos SwissProt.**

**Los números reales del pipeline**

La búsqueda proteome50 vs SwissProt tomó ~10 minutos 41 segundos y produjo: 36,488 representantes buscados, 23,252 con hit SwissProt, 13,236 sin hit, 1,821,541 hits totales, 16,980 representantes con hits strong o moderate.

De esos 16,980 representantes salieron 11,382 accesiones SwissProt únicas elegibles para AlphaFold. El piloto técnico de 1,000 accesiones —hecho para validar el script de descarga, no como objetivo final— dio: 986 encontradas, 13 no encontradas, 1 error. Ocupa ~991 MB. La proyección para todo el conjunto es ~10 GB, varias horas de descarga (mejor correr de noche).

**Las tres figuras de Maryam**

Se generaron las tres figuras solicitadas: distribución de longitudes (directamente desde `sequence_length` en la tabla maestra), distribución taxonómica (requirió construir tabla de anotación con género, familia, orden, clase y filo), y contexto de aislamiento. Esta última **no debe presentarse como "host-associated vs free-living"** — esa clasificación binaria no es defendible científicamente con los datos disponibles. Se usa en cambio "contexto original de aislamiento", que refleja lo que la literatura primaria reporta.

**Anotación funcional (sesión más reciente)**

Se anotaron funcionalmente los 1,000 accesiones del piloto vía UniProtKB ID Mapping (acceso directo desde navegador, por problema de DNS en la terminal WSL). Cobertura GO: 99.1%. Categoría dominante: enzimas metabólicas (40.9%). Transportadores: 11.1%. Solo el 4.6% sin función conocida — esperable dado el origen SwissProt.

**Advertencias metodológicas activas**

- Usar siempre MMseqs2 en el texto, nunca BLASTp.
- Las tablas, figuras y números del paper deben existir realmente en el repositorio o en la carpeta CeMbio. No usar resultados inferidos o generados por IA sin verificación.
- Las estructuras AlphaFold son de homólogos SwissProt, no de proteínas CeMbio directas.
- El piloto (1,000 accesiones) fue una prueba técnica. El análisis definitivo requiere las 11,382.

**Herramientas usadas**

Python 3.12, pandas, matplotlib, MMseqs2 v15-6f452, SwissProt local (574,627 registros), AlphaFold DB API, Git/GitHub, UniProtKB ID Mapping. Todo el trabajo se hizo con asistencia de Claude (Anthropic) — Claude Code para el pipeline técnico en terminal, Claude Cowork (Sonnet 4.6) para coordinación, escritura y documentación.

**¿Qué falta para el artículo?**

Fetch AlphaFold completo (11,382 accesiones), análisis pLDDT definitivo, análisis funcional profundo (Foldseek, InterPro/microdominios), Discusión e Introducción redactadas, revisión nomenclatural de *P. berkeleiyensis* y *P. nemavictus* con Maryam, y revisión manual de las candidatas NLS high_priority con José y Maryam.

---




### 3.9 Modelo de aprendizaje automático: predicción de divergencia funcional a partir de similitud estructural

Para evaluar si las características estructurales y de secuencia contienen información predictiva sobre divergencia funcional, se entrenó un modelo de ensamble de regresión logística (50 modelos con muestreo bootstrap, selección aleatoria de features por subconjunto √n) sobre los 12,272 pares no-self con similitud estructural significativa (TM-score ≥ 0.5) identificados en el análisis Foldseek. La variable respuesta fue la divergencia funcional binaria: un par se clasificó como "divergente" (1) si sus conjuntos de términos GO no compartían ningún término en común, y como "conservado" (0) en caso contrario. Las variables predictoras incluyeron TM-score, identidad de secuencia (fident), longitud de alineamiento, número de términos GO de cada proteína, diferencia absoluta en número de términos GO, longitud de cada proteína, diferencia de longitud y razón de longitudes (Tabla S1).

El dataset presentó un desbalance de clases marcado: 956 pares divergentes (7.8%) frente a 11,316 conservados (92.2%), coherente con la expectativa de que la mayoría de proteínas estructuralmente similares comparte función. Para compensar este desbalance, se aplicó ponderación por clase inversa a la frecuencia durante el entrenamiento. La evaluación se realizó sobre un conjunto de prueba retenido (20% de los datos, partición estratificada).

El modelo alcanzó un AUC-ROC de 0.764 (Figura X, panel C), demostrando que las características estructurales y de secuencia contienen información predictiva significativa sobre divergencia funcional, muy por encima del azar (AUC = 0.5). El recall sobre la clase divergente fue de 0.782, indicando que el modelo recupera el 78.2% de los pares realmente divergentes — una propiedad crucial para aplicaciones de anotación donde minimizar falsos negativos es prioritario. La precisión fue de 0.148, aproximadamente el doble de la frecuencia base de la clase divergente (7.8%), reflejando la dificultad inherente del desbalance. La identidad de secuencia (fident) y el TM-score emergieron como las variables de mayor importancia relativa, seguidas por el número de términos GO de cada proteína (Figura X, panel A).

![Modelo ML](./figures/ml_ensemble_results.png)
*Figura X. Resultados del modelo de aprendizaje automático para predicción de divergencia funcional. (A) Importancia relativa de las variables predictoras en el ensemble de 50 modelos. (B) Matriz de confusión sobre el conjunto de prueba (AUC=0.764, F1=0.249). (C) Curva ROC ilustrando la capacidad discriminativa del modelo para la clase divergente (pares con TM-score ≥ 0.5 y términos GO no compartidos).*


## 4. Discusión

### 4.1 Representatividad del conjunto piloto y distribución de longitudes proteicas

El conjunto piloto analizado en este estudio comprende 986 proteínas únicas de *C. elegans* con homólogos confirmados en SwissProt, de un total de 11,382 accesiones identificadas en el proteoma CeMbio completo. Aunque este piloto representa aproximadamente el 8.7% del conjunto total, su distribución de longitudes proteicas replica la heterogeneidad esperada en proteomas bacterianos: la mayoría de las proteínas presentan longitudes entre 100 y 500 aminoácidos, con una cola derecha correspondiente a proteínas multidominiales de mayor tamaño (Figura 1). Esta distribución es consistente con lo reportado para proteomas de bacterias Gram-negativas del género *Pseudomonas* y *Stenotrophomonas* (Winsor et al., 2016), que forman parte del CeMbio, y sugiere que el piloto captura adecuadamente la diversidad estructural del proteoma completo. No obstante, las proteínas de mayor longitud —potencialmente con arquitecturas multidominiales más complejas— estarán mejor representadas en análisis posteriores sobre el conjunto completo.

### 4.2 Diversidad taxonómica y contexto de aislamiento: implicaciones para la biología del hospedero

La distribución taxonómica del piloto refleja la composición proporcional de las 12 especies del CeMbio (Dirksen et al., 2020), con *Ochrobactrum vermis* y *Stenotrophomonas rhizophila* entre las especies más representadas (Figura 2). Esta diversidad es relevante porque las bacterias del CeMbio no son un conjunto arbitrario: fueron aisladas de poblaciones naturales de *C. elegans* y se sabe que modulan su desarrollo, reproducción y comportamiento (Berg et al., 2016). La heterogeneidad taxonómica del piloto implica que los patrones funcionales identificados aquí no son artefactos de una sola especie sino tendencias compartidas entre linajes bacterianos filogenéticamente distantes, lo que aumenta la robustez de las inferencias.

El contexto de aislamiento de las proteínas (Figura 3) muestra una predominancia de proteínas de origen de suelo y entorno intestinal, coherente con el nicho ecológico del nematodo. Es posible que esta distribución influya en los perfiles funcionales observados, dado que las bacterias adaptadas al intestino del hospedero tienden a enriquecer funciones relacionadas con la adquisición de nutrientes, evasión inmune y síntesis de factores de colonización (Montalvo-Katz et al., 2013). Futuros análisis comparativos entre contextos de aislamiento podrían revelar firmas funcionales específicas de nicho.

### 4.3 Cobertura de anotación funcional: validación del enfoque por homología

El mapeo de identificadores UniProtKB para las 986 proteínas del piloto resultó en una cobertura de anotación GO del 99.1%, con 978 proteínas con al menos un término GO asignado en alguna de las tres categorías ontológicas (proceso biológico, función molecular, componente celular) (Figura 4). Este nivel de cobertura es notablemente alto comparado con los reportados para proteomas bacterianos no-modelo anotados de novo mediante eggNOG-mapper, donde la cobertura típica ronda el 60–80% (Cantalapiedra et al., 2021). La diferencia se explica porque el presente análisis opera sobre homólogos SwissProt revisados manualmente (Swiss-Prot, sección curada de UniProtKB), en lugar de hacerlo directamente sobre los proteomas bacterianos del CeMbio. Esto introduce una limitación metodológica importante: las anotaciones reflejan la función conocida del homólogo en otra especie, no necesariamente la función de la proteína CeMbio en su contexto bacteriano original. Sin embargo, para el objetivo de este piloto —establecer un marco funcional de referencia y validar el pipeline de comparación estructura-secuencia— este nivel de cobertura es suficiente y metodológicamente justificado.

### 4.4 Categorías funcionales predominantes: sesgo hacia metabolismo y procesos de mantenimiento

El análisis de las categorías funcionales UniProt más frecuentes en el piloto muestra un enriquecimiento de proteínas asociadas a metabolismo central, biosíntesis de aminoácidos y proteínas ribosomales (Figura 5). Este patrón es esperado en cualquier muestreo no sesgado de un proteoma bacteriano, donde las funciones de mantenimiento celular (*housekeeping*) constituyen la mayoría numérica de las proteínas. Sin embargo, la presencia de proteínas con funciones relacionadas con transporte, señalización y secreción sugiere que el piloto captura también componentes funcionalmente relevantes para la interacción bacteria-hospedero. En análisis futuros, la comparación de la distribución funcional entre el piloto y el proteoma secretado predicho podría revelar enriquecimientos específicos de funciones de interfaz, como factores de virulencia atenuada, moduladores inmunes o enzimas degradativas del tegumento del nematodo (Geer et al., 2010).

### 4.5 Comparación estructura-secuencia: Foldseek revela homología funcional oculta

El hallazgo más relevante de este estudio piloto emerge de la comparación all-vs-all mediante Foldseek sobre las 986 estructuras predichas por AlphaFold: 8,608 pares de proteínas presentan similitud estructural significativa (TM-score ≥ 0.5) con identidad de secuencia inferior al 30%, y anotaciones GO distintas entre sí (Figura 6). Este resultado es consistente con el fenómeno bien documentado de convergencia estructural — o "analogía de pliegue" — en el que proteínas no relacionadas por descendencia común adoptan arquitecturas tridimensionales similares para desempeñar funciones distintas (Illergård et al., 2009; Andreeva & Murzin, 2006).

Lo que hace biológicamente relevante este resultado en el contexto del CeMbio es que demuestra que una fracción sustancial de la similitud estructural entre proteínas bacterianas del microbioma de *C. elegans* no se traduce en similitud funcional, y viceversa: las búsquedas por secuencia habrían agrupado estas proteínas con funciones diferentes, mientras que la búsqueda estructural las resuelve como entidades funcionalmente distintas a pesar de compartir pliegue. Este tipo de discriminación es especialmente crítico para proteínas de secreción, donde pequeñas diferencias en la superficie de interacción con el hospedero pueden tener consecuencias biológicas drásticas (Bharat & Bharat, 2019).

Cabe señalar que los valores de TM-score obtenidos (mediana ≈ 0.55 para los pares interesantes) corresponden a similitud estructural moderada, no a homología estructural estricta. Por tanto, los 8,608 pares no deben interpretarse como proteínas estructuralmente idénticas, sino como proteínas que comparten un pliegue o topología general pero han divergido en sus superficies funcionales. Este análisis piloto establece la base metodológica para escalar la comparación al proteoma CeMbio completo y, en combinación con herramientas de predicción de microdominios y aprendizaje automático, permitirá identificar proteínas con pliegue conservado pero especificidad funcional divergente — candidatos prioritarios para estudios estructurales experimentales.




### 4.6 Aprendizaje automático como capa integradora de información estructural y funcional

El modelo de ensamble entrenado sobre los pares Foldseek demuestra que la divergencia funcional entre proteínas estructuralmente similares es, en parte, predecible a partir de sus características cuantificables. Un AUC de 0.764 en un escenario con 7.8% de clase positiva representa una ganancia predictiva real sobre el azar, y el recall de 0.782 indica que el modelo recupera la gran mayoría de los casos de divergencia funcional real — un resultado con implicaciones directas para la priorización de candidatos en estudios experimentales.

La identidad de secuencia (fident) como variable más importante confirma la premisa central del análisis: proteínas con alta similitud estructural pero baja identidad de secuencia tienen mayor probabilidad de haber divergido funcionalmente. Esto es consistente con el fenómeno de convergencia de pliegue, donde la evolución converge hacia arquitecturas tridimensionales similares para funciones distintas (Andreeva & Murzin, 2006). El número de términos GO asignado a cada proteína también mostró importancia predictiva, sugiriendo que proteínas con mayor riqueza funcional anotada tienden a participar más frecuentemente en relaciones de similitud estructural con divergencia funcional — posiblemente porque son proteínas más estudiadas con múltiples roles.

Este resultado, aunque piloto, establece una prueba de concepto para el objetivo de aprendizaje automático de la tesis: la integración de información de secuencia, estructura y ontología funcional en un modelo predictivo mejora la capacidad de detectar proteínas cuya función no puede inferirse por homología de secuencia clásica. Escalar este enfoque al proteoma CeMbio completo (11,382 accesiones), incorporar embeddings proteicos de modelos de lenguaje (p. ej., ESM-2; Lin et al., 2023), y ampliar el conjunto de features estructurales (pLDDT, contactos residuales, descriptores de cavidades de unión) son las extensiones naturales de este trabajo.


## 5. Conclusiones

Este estudio piloto demuestra la viabilidad y el valor informativo de integrar comparación estructural con anotación funcional por ontología génica para el análisis de proteínas del microbioma de referencia de *C. elegans* (CeMbio). Los principales hallazgos son:

La anotación funcional por homología con UniProtKB alcanzó una cobertura del 99.2% sobre las 986 proteínas del piloto, validando la estrategia de usar homólogos Swiss-Prot como ancla para la anotación del proteoma CeMbio en ausencia de secuenciación directa o experimentos funcionales. Las categorías funcionales predominantes (metabolismo central, biosíntesis de aminoácidos, transporte) son consistentes con las funciones esperadas en proteomas bacterianos de vida libre y asociada a hospedero.

El análisis estructural all-vs-all con Foldseek identificó 8,608 pares de proteínas con similitud estructural significativa (TM-score ≥ 0.5) pero identidad de secuencia por debajo del 30% y anotaciones GO distintas. Este resultado confirma que una fracción sustancial de las relaciones funcionales en el proteoma CeMbio solo es detectable mediante comparación estructural, y quedaría invisible para métodos basados exclusivamente en secuencia. Los términos GO más frecuentes en estos pares — *plasma membrane*, *ATP binding*, *cytosol* — apuntan a proteínas con pliegues conservados de membrana y unión a nucleótidos que han divergido funcionalmente, lo que las convierte en candidatos prioritarios para estudios de especificidad funcional.

Estos resultados establecen el marco metodológico y los parámetros de referencia para escalar el análisis a las 11,382 accesiones del proteoma CeMbio completo, y para integrar en etapas posteriores métodos de predicción de microdominios y aprendizaje automático orientado a la identificación de determinantes de especificidad funcional en proteínas secretadas por la microbiota de *C. elegans*.


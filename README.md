nano README.md# celegans-struct-annotator

Local-first, modular software scaffold for the functional annotation of secreted bacterial proteins associated with *Caenorhabditis elegans*, with emphasis on structural similarity, microdomains, and reproducible downstream analysis.

## Project context

This repository is part of an undergraduate biology thesis focused on improving the functional annotation of secreted proteins from the microbiota associated with *C. elegans* by integrating sequence-based evidence, structural comparison, and later machine learning-based predictions.

The biological motivation is that many secreted proteins, especially short proteins or small functional regions, are difficult to annotate using sequence similarity alone. Structural comparison and integrative analysis can recover functional clues that traditional homology-based workflows often miss.

## Thesis title

**Functional annotation of microdomains in secreted proteins from the microbiota associated with *Caenorhabditis elegans* through the integration of sequence, structural, and machine learning information**

## Scientific motivation

Secreted bacterial proteins can modulate host processes such as immunity, development, longevity, and stress responses. However, functional annotation remains limited when proteins are highly divergent or contain short motifs, microdomains, or poorly characterized regions.

This project is motivated by three main ideas:

1. structure can remain informative even when sequence similarity is weak
2. microdomains and short functional regions are often underannotated
3. a reproducible software pipeline is needed to integrate local files, structural search, annotation comparison, and downstream analysis

## Software goal

The long-term goal of this repository is to become a clean, reproducible, extensible bioinformatics pipeline that can:

- load local protein metadata and sequence files
- validate and reuse local `.cif` structure files
- run Foldseek against configurable structural databases
- parse structural hits into clean tabular outputs
- compare structural evidence with functional annotation
- support downstream analysis and future API-assisted workflows

## Current status

This repository is currently in the **local MVP scaffold** stage.

At this point, the project is designed to support a local-first workflow and avoid hardcoded personal paths. It does **not** yet implement the full thesis pipeline.

### What is already scaffolded

- configuration loading from YAML
- local metadata loading from `.xlsx`, `.csv`, and `.tsv`
- local `.cif` inventory and validation
- safe Foldseek command preparation and validation
- Foldseek result parsing into normalized tables
- test suite for the MVP scaffold
- placeholders for future UniProt and AlphaFold integration

### What is intentionally deferred

- live UniProt API retrieval
- AlphaFold API / structure download integration
- large database automation
- annotation comparison against GO
- downstream statistical analysis and figure generation
- full machine learning integration

## Why local-first?

Previous attempts using large local structural resources such as AFDB50 / AlphaFold-UniProt50 became difficult to maintain because of:

- very large downloads
- unstable transfers
- corrupted temporary files
- confusing Foldseek database setup
- empty or misleading output files
- mixed WSL / Windows path issues

For that reason, this repository starts from a smaller, reproducible, local-first MVP and keeps large-scale database support as a later extension.

## Repository structure

```text
celegans-struct-annotator/
├── config/
├── data/
├── databases/
├── logs/
├── results/
├── scripts/
├── src/
├── structures/
├── tests/
├── README.md
└── pyproject.toml

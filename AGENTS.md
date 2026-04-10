# AGENTS

This repository is a local-first software scaffold for a biology thesis.

## Main rules
- Do not hardcode personal paths in source code.
- Prefer config-driven workflows.
- Reuse local CIF files before proposing downloads.
- Treat Domains as the current main workspace.
- Do not assume Gram is the current primary input.
- Do not rerun large Foldseek jobs unless explicitly requested.
- First validate parsing and CIF inventory using existing files.
- Keep the code modular, local-first, and reproducible.

## Current priority
Validate the real workflow with:
- local config loading
- existing CIF inventory
- existing non-empty Foldseek results
- clean tabular parsing

## Current preferred local inputs
- CIF directory: Domains/outputs/total_domains/afdb_cif
- Foldseek results: Domains/outputs/total_domains/results/test5_vs_swissprot.m8

## Avoid for now
- large reruns against SwissProt
- automatic API downloads
- premature ML integration
- mixing exploratory old folders into the new software core

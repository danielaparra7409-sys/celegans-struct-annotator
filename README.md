# celegans-struct-annotator

Local-first MVP scaffold for a bioinformatics thesis workflow focused on:

- loading local metadata from `.xlsx`, `.csv`, or `.tsv`
- validating local `.cif` files in a query directory
- running Foldseek against a configurable local target database
- parsing Foldseek results into a clean table

UniProt and AlphaFold integrations are placeholders only in this step.

## Scope

This scaffold only edits code, config, tests, and docs. It does not modify `data/`, `structures/`, `databases/`, `results/`, or `logs/`.

## Configuration

Copy [`config/config.example.yaml`](/home/daniela/code/celegans-struct-annotator/config/config.example.yaml) to a local config file and adjust the paths for your machine.

Supported metadata inputs:

- `.xlsx`
- `.csv`
- `.tsv`

## CLI

The MVP currently exposes exactly three entrypoints:

```bash
validate-local --config config/config.example.yaml
run-foldseek --config config/config.example.yaml
parse-results --config config/config.example.yaml
```

What each command does:

- `validate-local` loads metadata and inventories local `.cif` files from the configured query directory.
- `run-foldseek` validates the Foldseek setup, creates the configured temp directory if needed, and runs Foldseek with safe subprocess handling.
- `parse-results` reads the configured Foldseek TSV output and prints a normalized preview.

## Development

```bash
python -m pip install -e ".[dev]"
pytest
```

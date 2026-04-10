# Project State

## Current priority
Validate the software layer using the Domains workspace first.

## Main local source
Domains/

## Active CIF directory
Domains/outputs/total_domains/afdb_cif

## Active Foldseek results file for parser validation
Domains/outputs/total_domains/results/test5_vs_swissprot.m8

## Files currently considered unreliable or empty
- Domains/foldseek_run/results/all_vs_swissprot.m8
- Domains/foldseek_run/results/celegans_vs_swissprot.m8
- Domains/outputs/total_domains/foldseek_results/celegans_vs_swissprot.m8

## Current software goal
Validate the pipeline up to:
1. local config loading
2. CIF inventory and validation
3. parsing real Foldseek results
4. clean output table generation

## Not the goal right now
- rerun large SwissProt jobs
- full UniProt API integration
- AlphaFold API download automation
- full GO comparison
- ML integration

## Current local config
The project already uses `config/config.local.yaml`.

## Next step
Check CLI behavior with `config.local.yaml` and test `parse-results` using the non-empty `test5_vs_swissprot.m8` file.

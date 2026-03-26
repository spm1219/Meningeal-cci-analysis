# LIANA_root_scripts_for_CCI

This repository contains the root R scripts used to run LIANA-based cell–cell interaction (CCI) analysis from input count matrices and metadata files for the manuscript:

**Perilymphatic mast cell remodeling at the meningeal lymphatic interface in aged APP/PS2 mice**

## Files

- `01_run_liana_old_vs_young.R`  
  Reads `counts_OldvsYoung.csv`, `meta_OldvsYoung.csv`, and `Old_cluster_map.csv`, subsets the **Old** cells, and writes the aggregated LIANA result table as:
  - `liana_results_old.csv`

- `02_run_liana_nontg_vs_5xfad.R`  
  Reads `counts_NonTgvs5xFAD.csv`, `meta_NonTgvs5xFAD.csv`, and `AD_cluster_map.csv`, then runs LIANA separately for:
  - `liana_NonTg.csv`
  - `liana_5xFAD.csv`

- `Old_cluster_map.csv`  
  Cluster-to-cell-type annotation for the Old vs Young dataset.

- `AD_cluster_map.csv`  
  Cluster-to-cell-type annotation for the NonTg vs 5xFAD dataset.

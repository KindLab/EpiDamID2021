# EpiDamID2021

This repository contains key scripts and analyses accompanying the EpiDamID manuscript: [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.10.26.465688v1)

All raw data in the manuscript was processed using the workflow described in [Markodimitraki, Rang et al., 2020](https://www.nature.com/articles/s41596-020-0314-8), available at https://github.com/KindLab/scDamAndTools .

## Jupyter Notebooks
The following jupyter notebooks with analyses are included:
- LDA_classifier.applied_to_EB_data: Contains the analyses to train an LDA to predict membership to transcriptional clusters based on the single-cell DamID readout.

## Custom scripts
The following code is included:
- enrichment_over_regions.py: Custom code to compute average signal over regions of interest. Losely based on the DeepTools functionalities (https://deeptools.readthedocs.io/en/develop/).
- information_content.py: Custom code to compute the Information Content of DamID samples.

Additional code is available upon request.

# PAAD Genetic Landscape

Exploratory analysis of somatic mutations in pancreatic adenocarcinoma (TCGA-PAAD), focusing on the four most frequently altered driver genes: **KRAS, TP53, CDKN2A, SMAD4**. Built in R as a reproducible companion to my review on the genetic aspects of pancreatic ductal adenocarcinoma (PDAC).

## Data
- Source: cBioPortal, *Pancreatic Adenocarcinoma (TCGA, PanCancer Atlas)*
- File used: `data_mutations.txt` (not included here; download it from cBioPortal)
- 176 sequenced samples, 20,439 mutations (GRCh37)

## Methods
- R packages: `maftools`, `ggplot2`, `dplyr`
- Mutation summary, OncoPrint, lollipop plots, and pairwise somatic interaction tests (Fisher's exact test via `somaticInteractions`)

## Results
| Gene | Mutated samples | Frequency |
|---|---|---|
| KRAS | 117 | 66.5% |
| TP53 | 107 | 60.8% |
| SMAD4 | 37 | 21.0% |
| CDKN2A | 35 | 19.9% |

Figures are in `figures/`: mutation summary, OncoPrint, lollipop plots per gene, interaction plot, and mutation frequency bar chart.

## Limitations
- Frequencies reflect point mutations and indels only; copy-number alterations (e.g., CDKN2A homozygous deletions) are not included, so CDKN2A is likely underestimated.
- TCGA-PAAD contains samples with variable tumor purity, which can lower observed KRAS frequency compared with high-purity cohorts.

## How to run
1. Download `data_mutations.txt` from cBioPortal and place it next to the script.
2. Open `scripts/01_analysis.R` in RStudio and set the working directory to the folder containing the data.
3. Run the script. Figures are saved to `figures/`.

## Author
Alireza Hosseini, Laboratory Sciences, Shahid Beheshti University of Medical Sciences

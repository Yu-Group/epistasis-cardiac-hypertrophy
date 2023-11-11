# Code Repository for "Epistasis regulates genetic control of cardiac hypertrophy"

This repository contains code corresponding to:

Wang, Qianru, Tiffany M Tang, Nathan Youlton, Chad S Weldy, Ana M Kenney, Omer Ronen, J Weston Hughes, Elizabeth T Chin, Shirley C Sutton, Abhineet Agarwal, Xiao Li, Merle Behr, Karl Kumbier, Christine S Moravec, WH Wilson Tang, Kenneth B Margulies, Thomas P Cappola, Atul J Butte, Rima Arnaout, James B Brown, James R Priest, Victoria N Parikh, Bin Yu, Euan A Ashley. 2023. “[Epistasis regulates genetic control of cardiac hypertrophy](https://www.medrxiv.org/content/10.1101/2023.11.06.23297858v1).”

In summary, we developed an end-to-end pipeline to demonstrate the role of epistasis in the genetic control of cardiac hypertrophy. This pipeline consisted of four major phases: (1) the derivation of estimates of left ventricular mass via deep learning, (2) the computational prioritization of epistatic drivers based on UK Biobank via the low-signal signed iterative random forest (lo-siRF), (3) functional interpretation of the hypothesized epistatic genetic loci, and (4) experimental confirmation of epistasis in cardiac transcriptome and cardiac cellular morphology.

## Directory Structure

- **[01_lo-siRF/](./01_lo-siRF/)**: code to run low-signal signed iterative random forest (lo-siRF) given SNV and phenotype data as input
	- [dependencies.txt](./01_lo-siRF/dependencies.txt): package and dependency requirements to run lo-siRF pipeline
	- **[functions/](./01_lo-siRF/functions)**: contains helper functions to run lo-siRF analysis scripts
	- **[rmd/](./01_lo-siRF/rmd)**: contains files to reproduce [supplementary PCS R Markdown documentation](https://yu-group.github.io/epistasis-cardiac-hypertrophy/)
	- **[scripts/](./01_lo-siRF/scripts)**: contains scripts to run lo-siRF analysis (note: these scripts should be run in their numbered order via the corresponding \*.sh files)
		- [01_annovar.R](./01_lo-siRF/scripts/01_annovar.R): get SNV annotations (including the SNV to genomic locus mapping) using ANNOVAR
		- [02a_gwas_bolt_lmm.R](./01_lo-siRF/scripts/02a_gwas_bolt_lmm.R): run GWAS using BOLT-LMM as part of the dimension reduction stage in lo-siRF
		- [02b_gwas_plink.R](./01_lo-siRF/scripts/02b_gwas_plink.R): run GWAS using PLINK as part of the dimension reduction stage in lo-siRF
		- [03_prediction_models.R](./01_lo-siRF/scripts/03_prediction_models.R): fit common machine learning prediction models in R (e.g., siRF, RF, ridge, lasso) and assess prediction performance
		- [03_prediction_models.py](./01_lo-siRF/scripts/03_prediction_models.py): fit common machine learning prediction models in python (specifically SVM which is faster to train in python than R) and assess prediction performance
		- [04_sirf_local_stability.R](./01_lo-siRF/scripts/04_sirf_local_stability.R): fit siRF and the stability-driven importance measure to hypothesize and prioritize genomic loci and interactions between loci for follow-up experiments
		- [get_patient_characteristics.R](./01_lo-siRF/scripts/get_patient_characteristics.R): get summary statistics for the cohort under study (Supplementary Table 1)
		- [make_tables_and_figures.R](./01_lo-siRF/scripts/make_tables_and_figures.R): script to generate tables and figures in manuscript related to lo-siRF
- **[02_experimental_analysis/](./02_experimental_analysis/)**: contains functions and scripts to conduct statistical analysis of the microfluidics experimental data

## Citation

```r
@article{
	wang2023epistasis,
	author = {Qianru Wang and Tiffany M. Tang and Nathan Youlton and Chad S. Weldy and Ana M. Kenney and Omer Ronen and J. Weston Hughes and Elizabeth T. Chin and Shirley C. Sutton and Abhineet Agarwal and Xiao Li and Merle Behr and Karl Kumbier and Christine S. Moravec and W. H. Wilson Tang and Kenneth B. Margulies and Thomas P. Cappola and Atul J. Butte and Rima A. Arnaout and James B. Brown and James R. Priest and Victoria N. Parikh and Bin Yu and Euan A. Ashley},
	title = {Epistasis regulates genetic control of cardiac hypertrophy},
	year = {2023},
	doi = {10.1101/2023.11.06.23297858},
	publisher = {Cold Spring Harbor Laboratory Press},
	URL = {https://www.medrxiv.org/content/early/2023/11/09/2023.11.06.23297858},
	eprint = {https://www.medrxiv.org/content/early/2023/11/09/2023.11.06.23297858.full.pdf},
	journal = {medRxiv}
}
```

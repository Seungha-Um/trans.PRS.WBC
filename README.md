
# trans.PRS.WBC

Integrating genomic data with electronic health records (EHRs) requires multiple processing steps and specialized software for genomic data, particularly because of its complex structure and large scale. 
This repository provides workflows for computing a trans-ancestry PRS for WBC count (PRS-WBC) and linking the resulting scores to corresponding individual-level health records, along with downstream analyses.
We subsequently implement a phenome-wide association study (PheWAS) and tree-ensemble models to examine associations between genetic predisposition to WBC count and clinical phenotypes, as well as longitudinal WBC trajectories.

# Computing trans-ancestry PRS for WBC

### PRS-CSx

We applied PRS-CSx [2] to compute the trans-ancestry PRS-WBC. The PRS-CSx meethod is available at: https://github.com/getian107/PRScsx.

Implementation of PRS-CSx requires the following inputs:

1. Reference panel: We used the 1000 Genomes Project (1KG) reference panel. 
The reference panel is used (i) to identify common SNPs observed in both the target dataset and a well-genotyped reference dataset 
and (ii) to compute principal components (PCs) for the target data by projection onto the reference PC space. 
The publicly available 1KG data provide a standardized reference representing multiple ancestry groups.

2. GWAS summary statistics: GWAS summary statistics for each ancestry group were obtained from: https://www.mhi-humangenetics.org/resources/.
These summary statistics were generated from ancestry-specific meta-analyses combining multiple GWAS datasets [3].

3. LD reference panels: We used linkage disequilibrium (LD) reference panels constructed from the 1000 Genomes Project Phase 3 samples and provided: https://github.com/getian107/PRScsx

4. Target data: The target dataset consists of genotype data from individuals for whom PRS-WBC is to be calculated. In this study, we used genotype data from Phase III of the eMERGE Network.

Given these inputs, PRS-CSx generates ancestry-specific posterior SNP effect-size estimates by jointly modeling GWAS summary statistics from multiple ancestry groups. 
These posterior effect sizes account for the extent to which genetic effects are shared across ancestry groups while incorporating ancestry-specific LD patterns, GWAS sample sizes, 
and uncertainty in the original GWAS effect estimates. PRS-CSx uses Bayesian continuous shrinkage priors to shrink noisy SNP effect estimates while retaining variants with stronger evidence of association [2].

## Code in `inst/scripts`

### PRS computation

`PCA.R`: Performs PCA of the reference and target genotype data. The resulting PCs are used for post-hoc ancestry adjustment of PRS-WBC.

`PRS_CSx_computation_adjusted_for_PCs.R`: Computes PRS-WBC using the posterior SNP effect-size estimates obtained from PRS-CSx and performs post-hoc adjustment using genetic PCs.


`PRS_comparison_multiple_methods.R`: Compares the predictive performance of trans-ancestry PRS approaches with traditional ancestry-specific PRS approaches.

### PheWAS implementation

`create_wide_df_with_phecode_columns_for_PheWAS.R`: Creates a subject-level phenotype data frame containing phecode indicators derived from ICD-9 and ICD-10 diagnosis codes.

`PheWAS_implementation_and_figure_generation.R`: Performs the PheWAS using the generated phecode data and produces the corresponding figures and summary results.

### Longitudinal trajectory analysis using a tree-ensemble model

`regression_fit_BART_model.R`: Applies a Bayesian additive regression tree (BART) model with random intercept to estimate nonlinear longitudinal WBC trajectories while allowing for interactions among genetic and clinical risk factors.

`riBART_partial_dependency.R`: Generates the representative tree-split illustration presented in the manuscript.

### Summary statistics 

`checking ICD distribution_menopause_figure(Supp).R`: Examines the distribution of menopause-related ICD codes and generates the corresponding supplementary figure.

`WBC_median_computation.R`: Computes the median WBC count for each individual with repeated WBC measurements.

`Repeated_measures_summary_stats(Supp).R`: For the PheWAS implementation and PCA-based ancestry adjustment of PRS, we followed approaches described previously [1].


## Acknowledgments

The implementation of PheWAS and PCA-based ancestry adjustment of PRS in this repository was adapted from the code developed by Rosenthal et al. [1]. 
We gratefully acknowledge the authors for making their work and code available, which provided an important foundation for these analyses. 


For the longitudinal WBC trajectory analysis, we applied a random-intercept Bayesian additive regression tree (riBART) model to account for clustering across study sites [5] [6]. 
The `rbart_vi()` function in the `dbarts` R package is applicable. 


Alternatively, the `SoftBart` R package provides a flexible BART framework that can be extended to incorporate random intercepts, with the corresponding code provided by Linero et al. [4].
We gratefully acknowledge the authors for making their code publicly available, which provided an important foundation for these analyses.

## Reference 

[1] Rosenthal EA, Hsu L, Thomas M, Peters U, Kachulis C, Patterson K, Jarvik GP. Comparing ancestry standardization approaches for a transancestry colorectal cancer polygenic risk score. Genet Epidemiol. 2025;49(1). doi:10.1002/gepi.22590.

[2] Ruan Y, Lin YF, Feng YCA, et al. Improving polygenic prediction in ancestrally diverse populations. Nat Genet. 2022;54:573–580. doi:10.1038/s41588-022-01054-7.

[3] Chen MH, Raffield LM, Mousas A, et al. Trans-ethnic and ancestry-specific blood-cell genetics in 746,667 individuals from 5 global populations. Cell. 2020;182(5):1198–1213.e14. doi:10.1016/j.cell.2020.06.045.

[4] Linero AR, Yang Y. Bayesian regression tree ensembles that adapt to smoothness and sparsity. J R Stat Soc Series B Stat Methodol. 2018;80(5):1087–1110. doi:10.1111/rssb.12293.

[5] Tan YV, Roy J. Bayesian additive regression trees and the General BART model. Statistics and Its Interface. 2018;11(4):557–572.

[6] Tan YV, Flannagan CAC, Elliott MR. Predicting human-driving behavior to help driverless vehicles drive: random intercept Bayesian Additive Regression Trees. Statistics and Its Interface. 2018;11(4):557–572. doi: 10.4310/SII.2018.v11.n4.a1.

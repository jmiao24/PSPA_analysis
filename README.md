# POP-Inf_analysis

This repository contains the analysis code for POP-Inf paper: [Assumption-lean and data-adaptive post-prediction inference](https://arxiv.org/abs/2311.14220).

## Simulations

Simulation scripts are included in the folder `simulations`. Results and figures can be reproduced via the following three steps:

1. Run simulations for mean estimation and OLS coefficient estimation in `1-1_sim_mean.Rmd` and `1-2_sim_ols.Rmd`, respectively. These scripts will automatically create a subfolder `results` containing raw results.
2. Summarize simulation results into csv files: `2-1_sum_mean.Rmd` and `2-2_sum_ols.Rmd`.
3. Reproduce simulation figures in the manuscript: `3_make_figure.Rmd`.

## Real data analysis

Our GTEx real data analysis is based on the codes and data used in [Hypergraph factorization for multi-tissue gene expression imputation](https://www.nature.com/articles/s42256-023-00684-8) and [The impact of sex on gene expression across human tissues](https://www.science.org/doi/10.1126/science.aba3066?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed). We thank the authors for sharing the codes.

There are five steps for the pipeline:
1. Clean the GTEx data from data downloaded from [GTEx](https://www.gtexportal.org/home/): `1_clean_gtex.R`. 
2. Run the [HYFA](https://www.nature.com/articles/s42256-023-00684-8) to impute the gene expression: `2_run_hyfa.sh`.
3. Process the observed and imputed gene expression: `3_process_gtex.R`.
4. Run POP-Inf and compare it with alternative approaches to detect sex-biased genes with imputed gene expression: `4_popinf.R`.
5. Make the Figure3 in the POP-Inf paper: `5_make_figure.R`.

## Contact 

Please submit an issue or contact Jiacheng (jiacheng.miao@wisc.edu) or Xinran (xinran.miao@wisc.edu) for questions.

## Reference
```
@article{miao2023assumption,
  title={Assumption-lean and Data-adaptive Post-Prediction Inference},
  author={Miao, Jiacheng and Miao, Xinran and Wu, Yixuan and Zhao, Jiwei and Lu, Qiongshi},
  journal={arXiv preprint arXiv:2311.14220},
  year={2023}
}
```

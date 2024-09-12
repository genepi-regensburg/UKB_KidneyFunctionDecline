# UKB_KidneyFunctionDecline

This repository contains R scripts used for the analyses in the paper **Analyzing longitudinal trait trajectories using GWAS identifies genetic variants for kidney function decline**.

In particular:

- run_sevenApproaches: for running the seven approaches (*difference model*, *time model RI&RS*, *age model RI&RS*, *age model RI&RS uncorrelated*, *age model RI-only*, *BLUPs&LinReg*, and *age model RI&RS 350K*) on a subset of pre-selected genetic variants using the *lmer()* function from the *lme4* package.
- run_longGWAS: for running a longGWAS with the *age model RI&RS 350K* using the *glmmkin()* function from the *GMMAT* package and the *glmm.gei()* function from the *MAGEE* package.
- simulations: for running the simulations underlying the results in Table 2.
- simulations_helpers: helper functions for simulations.

The summary statistics of seven approaches on the 595 variants with known cross-sectional association with eGFR are provided in Supplementary Data 1.
The summary statistics of the longGWAS will be made available for download from http://www.genepi-regensburg.de/gwas-summary-statistics.

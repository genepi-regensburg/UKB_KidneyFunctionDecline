# UKB_KidneyFunctionDecline

This repository contains R scripts used for the analyses in the paper "Novel genetic variants identified for kidney function decline and insights into genome-wide association analyses using longitudinal trait trajectories".

In particular:
  - One R script for the seven approaches (difference model, time model RI&RS, age model RI&RS, age model RI&RS uncorrelated, age model RI-only, BLUPs&LinReg, and age model RI&RS 350K) using the lmer() function from the lme4 package.
  - One R script for the longGWAS with the age model RI&RS 350K using the glmmkin() function from the GMMAT package and the glmm.gei() function from the MAGEE package.

The summary statistics of seven approaches on the 595 variants with known cross-sectional association with eGFR are provided in Supplementary Table 5.
The summary statistics of the longGWAS will be made available for download from http://www.genepi-regensburg.de/gwas-summary-statistics.

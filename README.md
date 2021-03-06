# Meta-analysis-R-codes
This repository contains codes that can be easily adopted on different meta-analyses. 
* composite w rel: R function that produces a new dataset with composite correlation and reliabilities calculated based on multiple effect sizes and multiple reliability coefficients reported in the same sample for the same pair of relation. 
* CorrMatrix: R function that organizes Pearson correlation results (based on the corr.test function in Psych package) into an APA-format matrix.
* MetaCalculation: R function that enables (1) doing sample-size weighted meta-analysis based on the metafor package, (2) calculating the 95%CI for tau2, and (3) organizing all meta-analytic results into table format.
* psychmeta: R function that enables (1) doing psychometric meta-analysis weighted by sample size and corrected for measurement errors in both x and y; and (3) organizing all meta-analytic results into table format. 
* PubBias: R function that calculates and organizes test results for publication biases based on: (1) Egger's test; (2) PET-PEESE method; (3) Cumulative meta-analysis based on precision; (4) fixed-random trim and fill acompanied with contour-enhanced funnel plot. 
* resGroup: R function that calculates and organizes sub-group moderator analyses results (based on metafor package) into table format. 
* resWLS: R function that calculates and organizes WLS meta-regression results (based on metafor package) into table format (making it easier to export into csv or excel files).
* SubGroup power_post: R function for calculating post hoc power for fixed- and mixed-effects sub-group analyses for moderator in meta-analysis.
* SubGroup power_prior: R function for calculating a priori power for fixed- and mixed-effects subgroup analyses for moderator in meta-analysis. 
* WLS power_post: R function for calculating post hoc power for fixed- and mixed-effects meta-regression for moderators in meta-analysis. 

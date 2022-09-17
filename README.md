# Bayestrat Introduction
Bayestrat is a Bayesian genetic association test tool which identifies single-nucleotide polymorphisms (SNPs) associated with an outcome variable while using principle components (PCs) to account for population stratification. The current implementation supports continuous and discrete covariates and continuous only outcome variable. The package consists of three functions below. 
1. The bayestrat function implements the procedures for conducting SNP association analysis with population stratification correction. 
2. The pcCluster function is designed for data clustering analysis based on PC scores and provides a 3-dimensional visualization tool for the clustering results. 
3. The bayestratSummary function summarizes the output from bayestrat.

Bayestrat supports the following flexible model choices and advance computation.
1. The choices of Laplace (recommanded) and normal priors on the PC effect parameters;
2. Covergence diagnostics;
3. Parallel computation with multiple cores;
4. The package interfaces with C language and may call external software [PLINK](https://www.cog-genomics.org/plink/) if needed for large scale PC computation.

More detailed information annd intructions on using the package are provided in the package vignette.

# Method Introduction
Bayestrat supports a large number of PCs to fully correct for the underlying population structures, given that the sufficient number of PCs is unknown a priori. Utilizing the shrinkage Laplace priors, Bayestrat is able to filter out the irrelevant information, and select the highly confounded PCs. Bayestrat can be viewed as a compromise between the principle component regression of adding a few top PCs and the Bayesian version of linear mixed model (LMM). By shrinking the effects of irrelevant information, Bayestrat is able to achieve low type I error and high power.

Bayestrat first computes PCs based on a null data set (if PC data are not provided by users), then conducts association analysis for each SNP on the testing data set. Bayestrat provides inferences for effect sizes of SNPs, PCs, and other covariates based on the Markov Chain Monte Carlo (MCMC) samples from the posterior distribution.

More detailed information annd intructions on using the package are provided in the package vignette.

# Downloading
The Bayestrat package can be downloaded from Github using R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html).
```{r,warning=FALSE,message=FALSE,results=F}
library(devtools)
install_github("Zilu-Liu/Bayestrat@main")
```

If the dependencies are not automatically downloaded (especially for Windows users), please manually download these R packages as dependencies: [flashpcaR](https://github.com/gabraham/flashpca), data.table, rgl, coda.


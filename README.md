# Bayestrat
Bayesian genetic association tests for single-nucleotide polymorphisms (SNPs) while using principal components (PCs) to account for population stratification.
# Downloading
The Bayestrat package can be downloaded from Github using R package devtools ([link here](https://cran.r-project.org/web/packages/devtools/index.html)).
```{r,warning=FALSE,message=FALSE,results=F}
library(devtools)
install_github("Zilu-Liu/Bayestrat@main")
```

If the dependencies are not automatically downloaded (especially for Windows users), please manually download these R packages as dependencies: flashpcaR ([link here](https://github.com/gabraham/flashpca)), data.table, rgl, coda.


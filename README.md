# TissueDrug

Abstract
--------
Research in oncology traditionally focuses on specific tissue type from which the cancer develops. However, advances in high-throughput molecular profiling technologies have enabled the comprehensive characterization of molecular aberrations in multiple cancer types. It was hoped that these large-scale datasets would provide the foundation for a paradigm shift in oncology which would see tumors being classified by their molecular profiles rather than tissue types, but tumors with similar genomic aberrations may respond differently to targeted therapies depending on their tissue of origin. There is therefore a need to reassess the potential association between pharmacological response and tissue of origin for therapeutic drugs, and to test how these associations translate from preclinical to clinical settings.

In this paper, we investigate the tissue specificity of drug sensitivities in large-scale pharmacological studies and compare these associations to those found in drug clinical indications. Our meta-analysis of the four largest in vitro drug screening datasets indicates that tissue of origin is strongly associated with drug response. We identify novel tissue-drug associations, which may present new avenues for drug repurposing. One caveat is that the vast majority of the significant associations found in preclinical settings do not concur with clinical indications. Accordingly, our results call for more testing to find the root cause of the discrepancies between preclinical and clinical observations.


Citation
--------

Tissue specificity of in vitro drug sensitivity. Fupan Yao, Seyed Ali Madani Tonekaboni, Zhaleh Safikhani, Petr Smirnov, Nehme El-Hachem, Mark Freeman, Venkata Satya Kumar Manem, Benjamin Haibe-Kains. BioRxiv 2016, doi: http://dx.doi.org/10.1101/085357


Full Reproducibility of the Analysis Results
--------------------------------------------

We describe below how to fully reproduce the figures and tables reported in the paper

1.  Set up the software environment

2.  Run the R scripts

3.  Generate figures

Set up the software environment
-------------------------------

We developed and tested our analysis pipeline using R running on linux and Mac OSX platforms. The following is a copy of `sessionInfo()` from the development environment in R

```
 R version 3.3.0 (2016-05-03), x86_64-pc-linux-gnu

Base packages: base, datasets, graphics, grDevices, methods, parallel, stats, utils

Other packages: AIMS 1.4.0, Biobase 2.34.0, BiocGenerics 0.20.0, biomaRt 2.28.0, cluster 2.0.5, data.table 1.10.4, 
devtools 1.12.0, e1071 1.6-8, Formula 1.2-1, gdata 2.17.0, genefu 2.6.0, ggplot2 2.2.1, gplots 3.0.1, Hmisc 4.0-2, 
iC10 1.1.3, iC10TrainingData 1.0.1, igraph 1.0.1, lattice 0.20-34, limma 3.30.11, MASS 7.3-45, Matrix 1.2-8, 
mclust 5.2.2, minpack.lm 1.2-1, mRMRe 2.0.5, pamr 1.55, PharmacoGx 1.5.1, piano 1.14.5, preprocessCore 1.34.0, 
prodlim 1.5.9, qpcR 1.4-0, rgl 0.97.0, robustbase 0.92-7,
snow 0.4-2, snowfall 1.84-6.1, survcomp 1.24.0, survival 2.40-1, WriteXLS 4.0.0, xtable 1.8-2

Loaded via a namespace (and not attached): acepack 1.4.1, amap 0.8-14, AnnotationDbi 1.34.3, assertthat 0.1, 
backports 1.0.5, base64enc 0.1-3, BiocParallel 1.8.1, bitops 1.0-6, bootstrap 2015.2, caTools 1.17.1, 
celestial 1.3, checkmate 1.8.2, class 7.3-14, colorspace 1.3-2, DBI 0.5-1, DEoptimR 1.0-8, digest 0.6.12, 
downloader 0.4, fastmatch 1.1-0, fgsea 1.0.2, foreign 0.8-67, grid 3.3.0, gridExtra 2.2.1, gtable 0.2.0, 
gtools 3.5.0, htmlTable 1.9, htmltools 0.3.5, htmlwidgets 0.8, httpuv 1.3.3, IRanges 2.6.0, jsonlite 1.2, 
KernSmooth 2.23-15, knitr 1.15.1, latticeExtra 0.6-28, lava 1.4.7, lazyeval 0.2.0, lsa 0.73.1, magicaxis 2.0.0, 
magrittr 1.5, mapproj 1.2-4, maps 3.1.1, marray 1.52.0, memoise 1.0.0, mime 0.5, munsell 0.4.3, nnet 7.3-12, 
plotrix 3.6-4, plyr 1.8.4, R6 2.2.0, RANN 2.5, RColorBrewer 1.1-2, Rcpp 0.12.9, RCurl 1.95-4.8, relations 0.6-6, 
reshape2 1.4.2, rmeta 2.16, rpart 4.1-10, RSQLite 1.1-2, S4Vectors 0.10.1, scales 0.4.1, sets 1.0-16, shiny 1.0.0, 
slam 0.1-35, sm 2.2-5.4, SnowballC 0.5.1, splines 3.3.0, stats4 3.3.0, stringi 1.1.2, stringr 1.2.0, 
SuppDists 1.1-9.4, survivalROC 1.0.3, tibble 1.2, tools 3.3.0, withr 1.0.2, XML 3.98-1.5
```

All these packages are available on [CRAN](http://cran.r-project.org) or [Bioconductor](http://www.bioconductor.org)

All necessary packages have `library(<package>)` calls within the R scripts themselves, or the script assumes a previous script has been run and thus should have loaded nessesary packages. 

Running R Scripts in Repository and figure generation
-------------------------------
it is mandatory for all scripts to be run within one RStudio instance, as scripts will reference generated output variables from other script files

load CTRPv2, CCLE, GDSC, and gCSI using the `downloadPSet()` functions in PharmacoGx

run `GSEA_with_AUC.R` to get enrichment scores in one matrix variable called `combined1`

load the XML clinicaltrials.gov reader from `drugResultGetter.R` , and run `wordmine.R` to generate the `wordmine` variable referenced in diagram generation 

all figures are generated through their respective R file and should run independently with varying pdf generation in working directory set using `setWD()`

`figure7p1.R` and `figure7p2.R` generates the data for figure 7, which was then exported to circos for visualization

Figure 1 was manually created using Microsoft Word

supplementary file 3 is generated at the end of `GSEA_with_AUC.R`

supplementary figure 2 generated at the end of `supplementary file 2.R`


Docker
-------------------------------

Currently publicly hosted at `https://hub.docker.com/r/gosuzombie/tissuedrug/`. Running the command `docker pull gosuzombie/tissuedrug` will download the latest compiled image

Alternatively, the Dockerfile is also within the repo. Building the Docker image can be done using the command `docker build .` within the repository directory. 

Afterwards, the image id can be found using `docker images`

Finally, an interactive shell can be instantiated by running `docker run -it <docker image id>`



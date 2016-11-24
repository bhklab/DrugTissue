# TissueDrug

Abstract
--------
Research in oncology is traditionally focussed on a particular  tissue type from which the cancer under study develops. However, advances in high-throughput molecular profiling technologies enable the comprehensive characterization of molecular aberrations in multiple cancer types. It was hoped that this new data would provide the foundation for a paradigm shift in oncology which would see tumors being classified by their molecular profiles rather than tissue types. However, tumors with genomic aberrations may respond differently to targeted therapies depending on their tissue of origin. There is therefore a need to reassess the potential association between pharmacological response and tissue of origin for cytotoxic and targeted therapies, as well as how these associations translate from preclinical to clinical settings. In this paper, we investigate the tissue specificity of drug sensitivities in large-scale pharmacological studies and compared these associations to those found in clinical trials. Our meta-analysis of the four largest in vitro drug screening datasets indicates that tissue of origin is strongly predictive of drug response. Moreover, we identify novel tissue-drug associations, which may present exciting new avenues for drug repurposing. One caveat is that the vast majority of the significant associations found in preclinical settings do not concur with clinical observations. Accordingly, our results call for more testing to find the root cause of the discrepancies between preclinical and clinical observations.


Citation
--------

soon<sup>TM<sup>

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
R version 3.3.0 (2016-05-03)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets 
[7] methods   base  

other attached packages:
 [1] pvclust_2.0-0       vegan_2.4-0         lattice_0.20-33    
 [4] permute_0.9-0       survcomp_1.22.0     prodlim_1.5.7      
 [7] survival_2.39-4     piano_1.12.0        xlsx_0.5.7         
[10] xlsxjars_0.6.1      rJava_0.9-8         RColorBrewer_1.1-2 
[13] gplots_3.0.1        mgcv_1.8-12         nlme_3.1-128       
[16] ggplot2_2.1.0       reshape2_1.4.1      VennDiagram_1.6.17 
[19] futile.logger_1.4.1 PharmacoGx_1.1.6  

loaded via a namespace (and not attached):
 [1] gtools_3.5.0         lsa_0.73.1           slam_0.1-35         
 [4] sets_1.0-16          splines_3.3.0        colorspace_1.2-6    
 [7] SnowballC_0.5.1      marray_1.50.0        sm_2.2-5.4          
[10] magicaxis_1.9.4      BiocGenerics_0.18.0  lambda.r_1.1.7      
[13] plyr_1.8.4           lava_1.4.3           stringr_1.0.0       
[16] munsell_0.4.3        survivalROC_1.0.3    gtable_0.2.0        
[19] caTools_1.17.1       labeling_0.3         Biobase_2.32.0      
[22] parallel_3.3.0       Rcpp_0.12.5          KernSmooth_2.23-15  
[25] relations_0.6-6      scales_0.4.0         limma_3.28.11       
[28] gdata_2.17.0         rmeta_2.16           plotrix_3.6-2       
[31] bootstrap_2015.2     digest_0.6.9         stringi_1.1.1       
[34] SuppDists_1.1-9.2    tools_3.3.0          bitops_1.0-6        
[37] magrittr_1.5         cluster_2.0.4        futile.options_1.0.0
[40] MASS_7.3-45          Matrix_1.2-6         downloader_0.4      
[43] igraph_1.0.1 
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

Currently publicly hosted at `https://hub.docker.com/r/gosuzombie/tissuedrug/`

the Docker file is also within the repo. Building the Docker image can be done using the command `docker build .` within the repository directory. 

Afterwards, the image id can be found using `docker images`

Finally, an interactive shell can be instantiated by running `docker run -it <docker image id>`



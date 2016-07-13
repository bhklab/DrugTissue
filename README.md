# TissueDrug

Abstract
--------
Oncology is traditionally organized by disease sites, in other words, which tissue the cancer originates from. Advances in high-throughput molecular profiling technologies enable the comprehensive characterization of molecular aberrations in multiple cancer types. It was hoped that these new data would provide the foundation for paradigm shift in oncology where tumors would be better classified by their molecular profiles rather than tissue of origin. However, there is evidence that tumors with specific mutations or gene fusions respond differently to targeted therapies depending on their tissue type. There is therefore a need to assess the potential association between pharmacological response and tissue of origin for cytotoxic and targeted therapies and how these associations translate from preclinical to clinical settings. We investigate here the tissue specificity of drug sensitivities in large-pharmacological studies and compared these associations to clinical trials. Our meta-analysis of the four largest in vitro drug screening datasets indicates that tissue of origin is strongly predictive of drug response. Moreover, we identify novel tissue-drug associations, which could present an exciting new avenue for drug repurposing studies. However we found that the vast majority of the significant associations in preclinical settings do not concur with clinical observations. Our results call for more testing to find the root cause of the discrepancies between preclinical and clinical observations.


Citation
--------

you cant lol

Full Reproducibility of the Analysis Results
--------------------------------------------

We describe below how to fully reproduce the figures and tables reported in the paper

1.  Set up the software environment

2.  Run the R scripts

3.  Generate figures

Set up the software environment (needs to be updated)
-------------------------------

We developed and tested our analysis pipeline using R running on linux and Mac OSX platforms. To mimic our software environment the following R packages should  be installed:

* R version  3.0.1 (2013-05-16), x86_64-unknown-linux-gnu

* Base packages: base, datasets, graphics, grDevices, grid, methods, parallel,  splines, stats, utils

* Other packages: amap 0.8-7, Biobase 2.20.0,  BiocGenerics 0.6.0,  colorspace 1.2-2, GSA 1.03, MASS 7.3-26,  plotrix 3.4-7, prodlim 1.3.7,  survcomp 1.10.0,  survival 2.37-4, vcd 1.2-13,  WriteXLS 2.3.1,  xtable 1.7-1

* Loaded  via a namespace (and not attached): bootstrap 2012.04-0, epibasix  1.3, KernSmooth 2.23-10,  rmeta  2.16, SuppDists 1.1-9, survivalROC  1.0.3,  tools 3.0.1

All these packages are available on CRAN (http://cran.r-project.org) or Bioconductor (http://www.bioconductor.org)

Running R Scripts in Repository and figure generation
-------------------------------
it is mandatory for all scripts to be run within one RStudio instance, as scripts will reference generated output variables from other script files

load CTRPv2, CCLE, GDSC, and gCSI using the downloadPSet functions in PharmacoGx

run GSEA_with_AUC.R to get enrichment scores in one matrix, named combined1

load the XML clinicaltrials.gov reader from drugResultGetter.R , and run wordmine.R to generate a wordmine variable referenced in diagram generation 

all figures in their names figures.R files should run independently (Ctrl - A, run) , although any venn diagram generation doesnt seem to wipe the previous figure and just superimposes it

all figures write to the plot window in RStudio. A simple addition would be to add the call pdf("figure.pdf") and dev.off() around every figure call to write directly to pdf. 

to generate the circos plot export the edges matrix generated at the end of "cytoscape edges.R" and import into cytoscape.  

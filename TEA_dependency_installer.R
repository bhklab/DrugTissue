source("https://bioconductor.org/biocLite.R")

list.of.packages <- c("ggplot2", "XML", "mgcv", "reshape2", "grid", "gridExtra", "gplots","pheatmap","RColorBrewer", 
                      "VennDiagram", "gdata", "e1071", "xtable", "data.table", "snowfall", "utils",
                      "Hmisc", "WriteXLS", "rgl", "qpcR", "devtools", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.utstat.utoronto.ca/', dependencies = TRUE, clean = TRUE, Ncpus = 2, verbose = TRUE, quiet = TRUE)

biocLite(pkgs = c("Biobase", "survcomp", "piano", "PharmacoGx", "biomaRt", "preprocessCore", "genefu"), threads = 2)

dependencies <- c(list.of.packages, "Biobase", "survcomp", "piano", "PharmacoGx", "biomaRt", "preprocessCore", "genefu")


message("=== package install success status ===")

for(package in dependencies){
  message(paste(package, library(package, character.only = TRUE, quietly = TRUE, logical.return = TRUE, verbose = FALSE), sep = ": "))
}
source("https://bioconductor.org/biocLite.R")

list.of.CRAN.packages <- c("parallel", "ggplot2", "XML", "mgcv", "reshape2", "grid", "gridExtra", "gplots","pheatmap","RColorBrewer", 
                      "VennDiagram", "gdata", "e1071", "xtable", "data.table", "snowfall", "utils",
                      "Hmisc", "WriteXLS", "rgl", "qpcR", "devtools", "ggplot2")
new.CRAN.packages <- list.of.CRAN.packages[!(list.of.CRAN.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.CRAN.packages, repos='http://cran.utstat.utoronto.ca/', dependencies = TRUE, clean = TRUE, Ncpus = 2, verbose = TRUE, quiet = TRUE)


list.of.bioC.packages <- c("Biobase", "survcomp", "piano", "PharmacoGx", "biomaRt", "preprocessCore", "genefu")
new.bioC.packages <- list.of.bioC.packages[!(list.of.bioC.packages %in% installed.packages()[,"Package"])]
biocLite(pkgs = new.bioC.packages, Ncpus = 2, ask = FALSE)

dependencies <- c(list.of.CRAN.packages, list.of.bioC.packages)

message("=== package install success status ===")

for(package in dependencies){
  message(paste(package, library(package, character.only = TRUE, quietly = TRUE, logical.return = TRUE, verbose = FALSE), sep = ": "))
}

### Write.xls support for XLSX
installXLSXsupport()
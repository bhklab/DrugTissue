source("https://bioconductor.org/biocLite.R")

list.of.packages <- c("ggplot2", "XML", "mgcv", "reshape2", "grid", "gridExtra", "gplots", "pheatmap","RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(mgcv)
library(reshape2)
library(gridExtra)
library(grid)
library(gplots)
library(RColorBrewer)
library(openxlsx)
library(pheatmap)


if(!require(survcomp))
{
  biocLite("survcomp")  
  library(survcomp)
}

if(!require(piano))
{
  biocLite("piano")
  library(piano)
}
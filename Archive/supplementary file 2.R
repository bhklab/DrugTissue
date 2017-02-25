.libPaths("/rlibs/")
library(PharmacoGx)

if(!file.exists("./temp/combined1.RData"))
{
  source("./gsea_with_AUC.R")
} else {
  load("./temp/combined1.RData")
}

CCLE <- downloadPSet("CCLE_2013")
GDSC100 <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")
CTRPv2 <- downloadPSet("CTRPv2")



listofdrugs <- c(rownames(GDSC1000@drug), rownames(gCSI@drug), rownames(CCLE@drug), rownames(CTRPv2@drug))
finallist <- NULL

for(d in unique(listofdrugs))
{
  if(length(which(listofdrugs == d)) >= 3)
  {
    finallist <- c(finallist, d)
  }
}

for(d in finallist)
{
  out <- matrix(nrow = length(which(listofdrugs == d)) + 1, ncol = length(c("dataset", rownames(combined1))))
  colnames(out) <- c("dataset", rownames(combined1))
  counter <- 1
  for(ds in 1:length(ml))
  {
    if(d %in% colnames(ml[[ds]]))
    {
      for(a in 1:nrow(ml[[ds]]))
      {
        if(rownames(ml[[ds]])[a] %in% colnames(out))
        {
          out[counter, rownames(ml[[ds]])[a]] <- ml[[ds]][a, d]
          out[nrow(out), rownames(ml[[ds]])[a]] <- combined1[rownames(ml[[ds]])[a], d]
        }
        out[counter, "dataset"] <- names(ml)[ds]
      }
      counter <- counter + 1
    }
  }
  out[nrow(out), "dataset"] <- "combined"
  out[is.na(out)] <- ""
  write.xlsx(t(out), file = "./output/Supplementary_file_2.xlsx", sheetName = d, append = TRUE)
}
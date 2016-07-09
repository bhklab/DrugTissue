listofdrugs <- c(rownames(GDSC@drug), rownames(gCSI@drug), rownames(CCLE@drug), rownames(CTRPv2@drug))
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
  for(ds in ml)
  {
    if(d %in% colnames(ds))
    {
      for(a in 1:nrow(ds))
      {
        if(rownames(ds)[a] %in% colnames(out))
        {
          out[counter, rownames(ds)[a]] <- ds[a, d]
          out[nrow(out), rownames(ds)[a]] <- combined1[rownames(ds)[a], d]
        }
      }
      out[counter, "dataset"] <- names(ml)[counter]
      counter <- counter + 1
    }
  }
  out[nrow(out), "dataset"] <- "combined"
  write.xlsx(out, file = "suppfile1.xlsx", sheetName = d, append = TRUE)
}
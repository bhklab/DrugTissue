wordmine <- data.frame()

inputdruglist <- unique(c(rownames(GSK@drug), rownames(CCLE@drug), rownames(GDSC@drug), rownames(ctrpdrugs), breast, colorectal, prostate))
inputdruglist <- unique(toupper(gsub(badchars, "", inputdruglist)))

for(x in inputdruglist)
{
  if(! x %in% wordmine[1,])
  out <- getDrugTests(x)
  if(!is.na(out))
  {
    for(y in 1:nrow(out))
    {
      wordmine[x, out[y, "X2"]] <- out[y, "X1"]
    }
  }
}
save(wordmine, file = "./wordmine.RData")

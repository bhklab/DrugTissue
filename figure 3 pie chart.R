f <- NULL
for(x in rownames(combined1))
{
  for(y in colnames(combined1))
  {
    if(combined1[x,y] <= 0.05 && !is.na(combined1[x,y]))
    {
      f <- c(f, y)
    }
  }
}

length(unique(f))
#532
ncol(combined1)
#638

pie(c(523/638, 106/647), labels = c("specific: 523/638", "non-specific: 106/647"), cex = 2)
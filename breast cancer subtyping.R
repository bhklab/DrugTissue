library(plotrix)

breastsnosubtype <- combined["breast", ]
breastsubtype <- combined[grep("breast", rownames(combined)), ]

bsc <- breastsubtype

for(x in colnames(bsc))
{
  if(length(breastsnosubtype["breast", x]) > 0)
  {
    bsc["breast", x] <- breastsnosubtype["breast", x]
  }
  else
  {
    bsc["breast", x] <- NA
  }
  
}


ccc <- bsc

for(x in rownames(ccc))
{
  if(!(FALSE %in% is.na(ccc[x, ])))
  {
    ccc <- ccc[rownames(ccc) != x, ]
  }
}

for(x in colnames(ccc))
{
  if(!(FALSE %in% is.na(ccc[, x])))
  {
    ccc <- ccc[,colnames(ccc) != x ]
  }
}



my_palette <- colorRampPalette(c("#ece7f2", "#a6bddb", "#67a9cf"))



for(x in 1:nrow(ccc))
{
  for(y in 1:ncol(ccc))
  {
    if(!is.na(ccc[x,y]))
    {
      if(ccc[x,y] > 0.1)
      {
        ccc[x,y] <- "Not Significant"
      }
      else if(ccc[x,y] < 0.1 && ccc[x,y] > 0.05)
      {
        ccc[x,y] <- "< 0.1000"
      }
      else if(ccc[x,y] < 0.05 && ccc[x,y] > 0.01)
      {
        ccc[x,y] <- "< 0.05"
      }
      else if(ccc[x,y] < 0.01 && ccc[x,y] > 0.001)
      {
        ccc[x,y] <- "< 0.01"
      }
      else
      {
        ccc[x,y] <- "< 0.001"
      }
    }
    else if(is.nan(ccc[x,y]))
    {
      ccc[x,y] <- "Not Significant"
    }
  }
}
ccc <- as.matrix(ccc)

ccc <- as.data.frame(ccc)


cccm <- cccm[cccm$Var2 != "NVP-BEZ235",]

cccm <- ccc
rownames(cccm) <- toupper(rownames(cccm))
colnames(cccm) <- toupper(colnames(cccm))
cccm <- cccm[, intersect(toupper(gsub("-", ".",breast)), toupper(colnames(cccm)))]

cccm <- melt(as.matrix(cccm))
cccm$value <- as.factor(cccm$value)
cccm <- na.omit(cccm)
cccm$value <- as.factor(cccm$value)

ggplot(cccm, aes(cccm$Var2, cccm$Var1) ) + geom_tile(aes(fill = cccm$value), color = "black") + 
  scale_fill_manual(name = "p value", values = rev(brewer.pal(5,"Blues"))) + 
  labs(x = "drugs", y = "cellines") + 
  theme(axis.text.x = element_text(angle = 270, hjust = 1)) + coord_fixed(ratio = 1.5)

#pie charts 
dataset <- GDSC
size <- data.frame()
for(a in unique(dataset@cell$tissueid))
{
  if(!is.na(a))
  {
    size[a, "num"] <- length(grep(a, dataset@cell$tissueid))
  }
}

pie(size$num,labels = rownames(size), col=rainbow(nrow(size)), radius = , cex = 0.5)

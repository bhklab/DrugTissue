library(ggplot2)

path.input <- "Downloads"

breast <- xlsx::read.xlsx(paste(path.input, "12943_2015_312_MOESM2_ESM.xlsx", sep = "/"), sheetName = "Breast cancer")
colorectal <- xlsx::read.xlsx(paste(path.input, "12943_2015_312_MOESM2_ESM.xlsx", sep = "/"), sheetName = "Colorectal cancer")
prostate <- xlsx::read.xlsx(paste(path.input, "12943_2015_312_MOESM2_ESM.xlsx", sep = "/"), sheetName = "Prostate cancer")

breast <- breast[, "Drug.name"]
colorectal <- colorectal[, "Drug.name"]
prostate <- prostate[, "Drug.name"]

ccc <- combined
colnames(ccc) <- toupper(colnames(ccc))

for(a in breast)
{
  if(toupper(gsub("-", ".",a)) %in% colnames(ccc))
  {
    message(paste(a, "breast", sep = " "))
    print(ccc["breast", toupper(gsub("-", ".",a))])
  }
}

for(a in colorectal)
{
  if(toupper(gsub("-", ".",a)) %in% colnames(ccc))
  {
    message(paste(a, "large_intestine ", sep = " "))
    print(ccc["large_intestine", toupper(gsub("-", ".",a))])
  }
}

for(a in prostate)
{
  if(toupper(gsub("-", ".",a)) %in% colnames(ccc))
  {
    message(paste(a, "prostate", sep = " "))
    print(ccc["prostate", toupper(gsub("-", ".",a))])
  }
}

combinedcopy <- combined
ccc <- combinedcopy
ccc <- ccc[, h$colInd]

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
ccc <- as.numeric(ccc)
#ccc <- as.numeric(ccc)
heatmap.2(ccc,
          xlab = "drug",
          ylab =  "celllines",
          ccc,  # same data set for cell labels
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="both",     # only draw a row dendrogram
          lhei = c(1.5,10)) 

ccc <- as.data.frame(ccc)
rownames(ccc) <- toupper(rownames(ccc))
colnames(ccc) <- toupper(colnames(ccc))
ccc <- ccc[, intersect(toupper(gsub("-", ".", breast)), colnames(ccc))]


cccm <- melt(as.matrix(ccc))
cccm$value <- as.factor(cccm$value)
cccm <- na.omit(cccm)
cccm$value <- as.factor(cccm$value)

ggplot(cccm, aes(cccm$Var2, cccm$Var1) ) + geom_raster(aes(fill = cccm$value)) + 
  scale_fill_manual(values = rev(brewer.pal(5,"Blues"))) + 
  labs(x = "drugs", y = "cellines") + 
  theme(axis.text.x = element_text(angle = 270, hjust = 1))

cccm$Var1 <- toupper(cccm$Var1)
cccm$Var2 <- toupper(cccm$Var2)

ggplot(cccm, aes(cccm$Var1, cccm$Var2) ) + geom_tile(aes(fill = cccm$value), color = "black") + 
  scale_fill_manual(name = "p value", values = rev(brewer.pal(5,"Blues"))) + 
  labs(x = "cellines", y = "drugs") + 
  theme(axis.text.x = element_text(angle = 270, hjust = 1))  + coord_fixed(ratio = 1.5)
 


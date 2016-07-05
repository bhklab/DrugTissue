combinedcopy <- combined1[, c("Topotecan", "17-AAG", "Irinotecan", "paclitaxel", "Panobinostat")]

colnames(combinedcopy) <- toupper(gsub(badchars, "", colnames(combinedcopy)))
combinedcopy <- combinedcopy[grep("breast_", rownames(combinedcopy)) * -1, ]
ccc <- combinedcopy

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


ccc3 <- ccc
for(x in 1:nrow(ccc3))
{
  for(y in 1:ncol(ccc3))
  {
    if(!is.na(ccc3[x,y]))
    {
      if(ccc3[x,y] > 0.1)
      {
        ccc3[x,y] <- 0
      }
      else if(ccc3[x,y] < 0.1 && ccc3[x,y] > 0.05)
      {
        ccc3[x,y] <- 1
      }
      else if(ccc3[x,y] < 0.05 && ccc3[x,y] > 0.01)
      {
        ccc3[x,y] <- 2
      }
      else if(ccc3[x,y] < 0.01 && ccc3[x,y] > 0.001)
      {
        ccc3[x,y] <- 3
      }
      else
      {
        ccc3[x,y] <- 4
      }
    }
    else
    {
      ccc3[x,y] <- 0
    }
  }
}

my_palette <- colorRampPalette(c("#ece7f2", "#a6bddb", "#67a9cf"))

ccc3 <- as.matrix(ccc3)
h <- heatmap.2(ccc3,
               xlab = "drug",
               ylab =  "celllines",
               ccc3,  # same data set for cell labels
               density.info="histogram",  # turns off density plot inside color legend
               trace="none",         # turns off trace lines inside the heat map
               col=my_palette,       # use on color palette defined earlier 
               dendrogram="both",     # only draw a row dendrogram
               lhei = c(1.5,10))

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
        ccc[x,y] <- "< 0.10"
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

cccm <- melt(as.matrix(ccc))
cccm$value <- as.factor(cccm$value)
cccm <- na.omit(cccm)
cccm$value <- as.factor(cccm$value)

ggplot(cccm) + geom_tile(aes(fill = cccm$value, x = cccm$Var1, y = cccm$Var2), color = "black") + 
  scale_fill_manual(name = "p value", values = rev(brewer.pal(5,"Oranges")),limits = rev(c("Not Significant", "< 0.10", "< 0.05", "< 0.01", "< 0.001")), drop = FALSE) + 
  labs(x = "", y = "") + 
  scale_y_discrete(limits = colnames(ccc)[h$colInd]) +
  scale_x_discrete(limits = rownames(ccc)[h$rowInd]) +
  theme(axis.text.x = element_text(angle = 270, hjust = 1, size = 15), axis.text.y = element_text(size = 20))
#normalized by tissue amount
graph <- data.frame(matrix(nrow = 0,ncol = 3))
colnames(graph) <- c("dataset", "category", "value")

for(ds in 1:length(ml))
{
  imp <- ml[ds]
  imp <- imp[[1]]
  
  dscounter <- 0
  metacounter <- 0
  for(x in rownames(imp))
  {
    for(y in colnames(imp))
    {
      if(!is.na(imp[x,y]) && imp[x,y] <= 0.05)
      {
        dscounter <- dscounter + 1
        
        if(y %in% colnames(combined1))
        {
          if(!is.na(combined1[x,y]) && combined1[x,y] <= 0.05)
          {
            metacounter <- metacounter + 1
          }
        }
        
      }
    }
  }
  
  dscounter <- dscounter - metacounter
  
  if(names(ml)[ds] == "GDSC1000")
  {
    dscounter <- dscounter/sum(linenums$GDSC1000)
  }
  if(names(ml)[ds] == "CCLE")
  {
    dscounter <- dscounter/sum(linenums$CCLE)
  }
  if(names(ml)[ds] == "gSCI")
  {
    print(dscounter)
    print(sum(linenums$gCSI))
    dscounter <- dscounter/sum(linenums$gCSI)
  }
  if(names(ml)[ds] == "CTRP")
  {
    dscounter <- dscounter/sum(linenums$CTRPv2)
  }
  
  metacounter <- metacounter/sum(linenums[!is.na(linenums)])
  
  graph <- rbind(graph, data.frame(dataset=names(ml)[ds], category="metaanalysis", value=metacounter))
  graph <- rbind(graph, data.frame(dataset=names(ml)[ds], category="dataset", value=dscounter))
}

pdf("figure21celllines.pdf")
ggplot(graph, aes(dataset, value)) + 
  geom_bar(aes(fill = reorder(graph$category, rev(graph$value))), stat = "identity", width = 0.8) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20)) +
  guides(fill=guide_legend(title="")) +
  scale_x_discrete(limits = c("CCLE", "gSCI", "GDSC1000", "CTRP")) + 
  labs(x = NULL, y = "ratio of tissuedrug interactions to cell lines") +
  theme(legend.text=element_text(size=10)) +
  coord_fixed(ratio = 4, expand = TRUE)
dev.off()


#normalized by drugs
graph <- data.frame(matrix(nrow = 0,ncol = 3))
colnames(graph) <- c("dataset", "category", "value")

for(ds in 1:length(ml))
{
  imp <- ml[ds]
  imp <- imp[[1]]
  
  dscounter <- 0
  metacounter <- 0
  for(x in rownames(imp))
  {
    for(y in colnames(imp))
    {
      if(!is.na(imp[x,y]) && imp[x,y] <= 0.05)
      {
        dscounter <- dscounter + 1
        
        if(y %in% colnames(combined1))
        {
          if(!is.na(combined1[x,y]) && combined1[x,y] <= 0.05)
          {
            metacounter <- metacounter + 1
          }
        }
        
      }
    }
  }
  
  dscounter <- dscounter - metacounter
  dscounter <- dscounter/ncol(ml[ds])

  
  metacounter <- metacounter/ncol(combined1)
  
  graph <- rbind(graph, data.frame(dataset=names(ml)[ds], category="metaanalysis", value=metacounter))
  graph <- rbind(graph, data.frame(dataset=names(ml)[ds], category="dataset", value=dscounter))
}

pdf("figure21drugs.pdf")
ggplot(graph, aes(dataset, value)) + 
  geom_bar(aes(fill = reorder(graph$category, rev(graph$value))), stat = "identity", width = 0.8) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=20)) +
  guides(fill=guide_legend(title="")) +
  scale_x_discrete(limits = c("CCLE", "gSCI", "GDSC1000", "CTRP")) + 
  labs(x = NULL, y = "ratio of tissuedrug interactions to drugs") +
  theme(legend.text=element_text(size=10)) +
  coord_fixed(ratio = 4, expand = TRUE)
dev.off()




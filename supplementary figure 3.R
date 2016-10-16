graph <- data.frame(matrix(nrow = 0,ncol = 3))
colnames(graph) <- c("dataset", "category", "value")

for(ds in 1:length(ml))
{
  imp <- ml[ds]
  imp <- imp[[1]]
  
  dscounter <- 0
  metacounter <- 0
  dsicounter <- 0
  for(x in rownames(imp))
  {
    for(y in colnames(imp))
    {
      if(!is.na(imp[x,y]) && imp[x,y] < 0.05)
      {
        dscounter <- dscounter + 1
        
        if(y %in% colnames(combined1))
        {
          if(!is.na(combined1[x,y]) && combined1[x,y] < 0.05)
          {
            metacounter <- metacounter + 1
          }
        }
      }
      else if(!is.na(imp[x,y] && imp[x,y] > 0.05))
      {
        if(!is.na(combined1[x,y]) && combined1[x,y] < 0.05)
        {
          dsicounter <- dsicounter + 1
        }
      }
    }
  }
  
  dscounter <- dscounter - metacounter
  graph <- rbind(graph, data.frame(dataset=names(ml)[ds], category="metaanalysis", value=metacounter))
  graph <- rbind(graph, data.frame(dataset=names(ml)[ds], category="dataset", value=dscounter))
  graph <- rbind(graph, data.frame(dataset=names(ml)[ds], category="insignificant", value=dsicounter))
}

graph$dataset <- gsub("gSCI", "gCSI", graph$dataset)
graph <- rbind(graph, data.frame(dataset="metaanalysis", category="metaanalysis", value=length(combined1[!is.na(combined1) & combined1 < 0.05])))

pdf("Supplementary_Figure_3.pdf", width = 10, height = 10)
ggplot(graph, aes(dataset, value)) + 
  geom_bar(aes(fill = reorder(graph$category, rev(graph$value))), stat = "identity", width = 0.8) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=30)) +
  guides(fill=guide_legend(title="", keywidth=0.1, keyheight=0.5, default.unit="inch")) +
  scale_x_discrete(limits = c("gCSI", "CCLE", "CTRPv2", "GDSC1000", "metaanalysis")) + 
  #scale_y_log10() + 
  labs(x = NULL, y = "number of drug tissue interactions") +
  theme(legend.text=element_text(size=30))
dev.off()




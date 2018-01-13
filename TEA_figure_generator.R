### File Inputs
p_adj_false <- read.xlsx(file.path(GSEADir, "DrugTissueAssocs_All.xlsx"), sheet = 1)
p_adj_true <- read.xlsx(file.path(GSEADir, "DrugTissueAssocs_All.xlsx"), sheet = 2)

p_adj_true <- as.data.frame(p_adj_true)
p_adj_false <- as.data.frame(p_adj_false)

cut_off <- FDRcutoff

sig_adj_false <- p_adj_false[!is.na(p_adj_false$Combined_FDR) & as.numeric(p_adj_false$Combined_FDR) < cut_off, ]
sig_adj_true <- p_adj_true[!is.na(p_adj_true$Combined_FDR) & as.numeric(p_adj_true$Combined_FDR) < cut_off, ]

### Library Load
library(PharmacoGx)
library(VennDiagram)
library(ggplot2)
library(xlsx)

GDSC1000 <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")
CCLE <- downloadPSet("CCLE")
CTRPv2 <- downloadPSet("CTRPv2") 

if(!dir.exists("./figures")){
  dir.create("./figures")
}


### Cell Line Breakdown by Dataset (Figure 2A)

linenums <- data.frame()

t <- table(CCLE@cell$tissueid)
linenums[names(t), "CCLE"] <- t

t <- table(GDSC1000@cell$tissueid)
linenums[names(t), "GDSC1000"] <- t

t <- table(gCSI@cell$tissueid)
linenums[names(t), "gCSI"] <- t

t <- table(CTRPv2@cell$tissueid)
linenums[names(t), "CTRPv2"] <- t

linenums[is.na(linenums)] <- 0
linenums <- linenums[rownames(linenums) != "",]

linenums$median <- apply(linenums, MARGIN = 1, median)
linenums$sum <- apply(linenums, MARGIN = 1, sum)
linenums <- linenums[linenums$sum > 8,]
linenums <- linenums[order(linenums$median), ]

linenums <- linenums[, c(-5, -6)]
order <- rownames(linenums)

linenums_melt <- melt(as.matrix(linenums))

colnames(linenums_melt) <- c("tissue type", "dataset", "value")

linenums_melt$`tissue type` <- gsub("_", " ", linenums_melt$`tissue type`)
linenums_melt$`tissue type` <- as.factor(linenums_melt$`tissue type`)

pdf(file.path(figure_output, "figure2a.pdf"))
ggplot(linenums_melt, aes(x = linenums_melt$`tissue type`, y = linenums_melt$value, fill=linenums_melt$dataset)) + 
  geom_bar(stat = "identity", position="dodge") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.title=element_text(size=10)) +
  guides(fill=guide_legend(title="")) +
  scale_x_discrete(limits = gsub("_", " ", order)) + 
  labs(x = NULL, y = "number of cell lines") +
  theme(legend.text=element_text(size=10))
dev.off()

### Venn Diagrams for Overlap (Figure 2B)

drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="N-{3-Chloro-4-[(3-fluorobenzyl)oxy]phenyl}-6-[5-({[2-(methylsulfonyl)ethyl]amino}methyl)-2-furyl]-4-quinazolinamine", rep="lapatinib", fixed=TRUE)
drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="bleomycin A2", rep="Bleomycin", fixed=TRUE)
drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="cytarabine hydrochloride", rep="Cytarabine", fixed=TRUE)
drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="doxorubicin", rep="Doxorubicin", fixed=TRUE)
drugNames(CCLE) <- gsub(drugNames(CCLE), pat = "AZD6244", rep = "selumetinib", fixed = TRUE)

ccled <- colnames(ml$CCLE)
gdscd <- colnames(ml$GDSC1000)
ctrpd <- colnames(ml$CTRPv2)
gcsid <- colnames(ml$gCSI)

pdf("./figures/figure2drugs.pdf", width = 100, height = 10, paper = "USr")
#figure 2 drugs
draw.quad.venn(
  area1 = length(ccled),
  area2 = length(gdscd),
  area3 = length(ctrpd),
  area4 = length(gcsid),
  n12 = length(intersect(ccled, gdscd)),
  n13 = length(intersect(ccled, ctrpd)),
  n14 = length(intersect(ccled, gcsid)),
  n23 = length(intersect(gdscd, ctrpd)),
  n24 = length(intersect(gdscd, gcsid)),
  n34 = length(intersect(ctrpd, gcsid)),
  n123 = length(intersect(intersect(ccled, gdscd), ctrpd)),
  n124 = length(intersect(intersect(ccled, gdscd), gcsid)),
  n134 = length(intersect(intersect(ccled, ctrpd), gcsid)),
  n234 = length(intersect(intersect(gdscd, ctrpd), gcsid)), 
  n1234 = length(intersect(intersect(ccled, gdscd), intersect(ctrpd, gcsid))), 
  category = c("CCLE", "GDSC1000", "CTRPv2", "gCSI"),
  cat.cex = 2,
  cex = 2,
  fill = c("#F8766D", "#B79F00", "#00BFC4", "#C77CFF"),
  fontfamily = rep("Helvetica", 15),
  cat.fontfamily = rep("Helvetica", 4)
)
dev.off()

cclecell <- CCLE@cell
gdsc1000cell <- GDSC1000@cell
ctrpv2cell <- CTRPv2@cell
gcsicell <- gCSI@cell

pdf("./figures/figure2cells.pdf", width = 100, height = 10, paper = "USr")
#figure 2 cells
draw.quad.venn(
  area1 = nrow(cclecell),
  area2 = nrow(gdsc1000cell),
  area3 = nrow(ctrpv2cell),
  area4 = nrow(gcsicell),
  n12 = length(intersect(rownames(cclecell), rownames(gdsc1000cell))),
  n13 = length(intersect(rownames(cclecell), rownames(ctrpv2cell))),
  n14 = length(intersect(rownames(cclecell), rownames(gcsicell))),
  n23 = length(intersect(rownames(gdsc1000cell), rownames(ctrpv2cell))),
  n24 = length(intersect(rownames(gdsc1000cell), rownames(gcsicell))),
  n34 = length(intersect(rownames(ctrpv2cell), rownames(gcsicell))),
  n123 = length(intersect(intersect(rownames(cclecell), rownames(gdsc1000cell)), rownames(ctrpv2cell))),
  n124 = length(intersect(intersect(rownames(cclecell), rownames(gdsc1000cell)), rownames(gcsicell))),
  n134 = length(intersect(intersect(rownames(cclecell), rownames(ctrpv2cell)), rownames(gcsicell))),
  n234 = length(intersect(intersect(rownames(gdsc1000cell), rownames(ctrpv2cell)), rownames(gcsicell))), 
  n1234 = length(intersect(intersect(rownames(cclecell), rownames(gdsc1000cell)), intersect(rownames(ctrpv2cell), rownames(gcsicell)))), 
  category = c("CCLE", "GDSC1000", "CTRPv2", "gCSI"),
  cat.cex = 2, 
  cex = 2,
  fill = c("#F8766D", "#B79F00", "#00BFC4", "#C77CFF"),
  fontfamily = rep("Helvetica", 15),
  cat.fontfamily = rep("Helvetica", 4)
)
dev.off()

gdsc1000cell$tissueid <- gsub("esophagus", "oesophagus", gdsc1000cell$tissueid, fixed = TRUE)
gdsc1000cell <- gdsc1000cell[gdsc1000cell$tissueid != "misc", ]
gdsc1000cell <- gdsc1000cell[gdsc1000cell$tissueid != "other", ]
ctrpv2cell <- ctrpv2cell[ctrpv2cell$tissueid != "", ]
gcsicell <- gcsicell[gcsicell$tissueid != "", ]

pdf("./figures/figure2tissue.pdf", width = 100, height = 10, paper = "USr")
#figure 2 tissue types
draw.quad.venn(
  area1 = length(unique(na.omit(cclecell$tissueid))),
  area2 = length(unique(na.omit(na.omit(gdsc1000cell$tissueid)))),
  area3 = length(unique(na.omit(ctrpv2cell$tissueid[ctrpv2cell$tissueid != ""]))),
  area4 = length(unique(na.omit(gcsicell$tissueid[gcsicell$tissueid != ""]))),
  n12 = length(intersect(na.omit(cclecell$tissueid), na.omit(gdsc1000cell$tissueid))),
  n13 = length(intersect(na.omit(cclecell$tissueid), ctrpv2cell$tissueid[ctrpv2cell$tissueid != ""])),
  n14 = length(intersect(na.omit(cclecell$tissueid), gcsicell$tissueid[gcsicell$tissueid != ""])),
  n23 = length(intersect(na.omit(gdsc1000cell$tissueid), ctrpv2cell$tissueid[ctrpv2cell$tissueid != ""])),
  n24 = length(intersect(na.omit(gdsc1000cell$tissueid), gcsicell$tissueid[gcsicell$tissueid != ""])),
  n34 = length(intersect(ctrpv2cell$tissueid[ctrpv2cell$tissueid != ""], gcsicell$tissueid[gcsicell$tissueid != ""])),
  n123 = length(intersect(intersect(na.omit(cclecell$tissueid), na.omit(gdsc1000cell$tissueid)), ctrpv2cell$tissueid[ctrpv2cell$tissueid != ""])),
  n124 = length(intersect(intersect(na.omit(cclecell$tissueid), na.omit(gdsc1000cell$tissueid)), gcsicell$tissueid[gcsicell$tissueid != ""])),
  n134 = length(intersect(intersect(na.omit(cclecell$tissueid), ctrpv2cell$tissueid[ctrpv2cell$tissueid != ""]), gcsicell$tissueid[gcsicell$tissueid != ""])),
  n234 = length(intersect(intersect(na.omit(gdsc1000cell$tissueid), ctrpv2cell$tissueid[ctrpv2cell$tissueid != ""]), gcsicell$tissueid[gcsicell$tissueid != ""])), 
  n1234 = length(intersect(intersect(na.omit(cclecell$tissueid), na.omit(gdsc1000cell$tissueid)), intersect(ctrpv2cell$tissueid[ctrpv2cell$tissueid != ""], gcsicell$tissueid[gcsicell$tissueid != ""]))), 
  category = c("CCLE", "GDSC1000", "CTRPv2", "gCSI"),
  cex = 2,
  cat.cex = 2,
  fill = c("#F8766D", "#B79F00", "#00BFC4", "#C77CFF"),
  fontfamily = rep("Helvetica", 15),
  cat.fontfamily = rep("Helvetica", 4)
)
dev.off()


### Interaction Bar Plot by Tissue Type (Figure 3)

graph <- data.frame(matrix(nrow = 0,ncol = 3))
colnames(graph) <- c("tissue", "drug", "value")

graph <- sig_adj_true[, c("Tissue", "Drug")]
graph$value <- 2
colnames(graph) <- c("tissue", "drug", "value")

for(entry in 1:nrow(sig_adj_false)){
  if(nrow(graph[graph$tissue == sig_adj_false$Tissue[entry] & graph$drug == sig_adj_false$Drug[entry],]) == 0)
  {
    graph <- rbind(graph, data.frame(tissue=sig_adj_false$Tissue[entry], drug=sig_adj_false$Drug[entry], value=1))
  } else {
    graph[graph$tissue == sig_adj_false$Tissue[entry] & graph$drug == sig_adj_false$Drug[entry], "value"] <- 3
  }
}

graph_real <- data.frame(matrix(nrow = 0,ncol = 3))
colnames(graph_real) <- c("tissue", "category", "value")

for(tissue in unique(graph$tissue)){
  graph_real <- rbind(graph_real, data.frame(tissue=tissue, category="original", value=nrow(graph[graph$tissue == tissue & graph$value == 1,])))
  graph_real <- rbind(graph_real, data.frame(tissue=tissue, category="overlap", value=nrow(graph[graph$tissue == tissue & graph$value == 3,])))
  graph_real <- rbind(graph_real, data.frame(tissue=tissue, category="adjusted", value=nrow(graph[graph$tissue == tissue & graph$value == 2,])))
}

graph_real$tissue <- gsub("_", " ", graph_real$tissue)

graph_real <- na.omit(graph_real)

pdf("./figures/figure3.pdf", width = 10)
ggplot(graph_real, aes(tissue, value)) + 
  geom_bar(aes(fill = graph_real$category), stat = "identity", width = 0.8) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=30)) +
  guides(fill=guide_legend(title="", keywidth=0.1, keyheight=0.5, default.unit="inch")) +
  scale_x_discrete(limits = gsub("_", " ", names(sort(table(graph$tissue))))) + 
  #scale_y_log10() + 
  labs(x = "tissue type", y = "number of interactions") +
  theme(legend.text=element_text(size=30))
dev.off()

### Interaction Bar Plot by Dataset (Figure 4)
list_of_sig_interactions <- unique(rownames(p_adj_false), rownames(p_adj_true))

graph <- data.frame(matrix(nrow = length(list_of_sig_interactions),ncol = 5))
rownames(graph) <- list_of_sig_interactions
colnames(graph) <- c("FDR_CCLE","FDR_gCSI", "FDR_CTRPv2", "FDR_GDSC1000", "Combined_FDR")
graph[is.na(graph)] <- 0

for(intact in rownames(graph)){
  if(intact %in% rownames(p_adj_true)){
    graph[intact, ] <- p_adj_true[intact, c("FDR_CCLE","FDR_gCSI", "FDR_CTRPv2", "FDR_GDSC1000", "Combined_FDR")]
    graph[intact,] <- ifelse(as.numeric(graph[intact,]) < cut_off & graph[intact,] == 0, 1,0)
  }
  if(intact %in% rownames(p_adj_false)){
    for(ps in colnames(graph)){
      if(!is.na(p_adj_false[intact, ps])){
        graph[intact, ps] <- ifelse(as.numeric(p_adj_false[intact, ps]) < cut_off && graph[intact, ps] == 0, 1, 0)
      }
    }
  }
}

graph_real <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(graph_real) <- c("dataset", "category", "value")

for(dataset in c("FDR_CCLE","FDR_gCSI", "FDR_CTRPv2", "FDR_GDSC1000")){
  both <- 0
  in_meta <- 0
  in_dataset <- 0
  
  for(r in rownames(graph)){
    
    if(is.na(graph[r, dataset])){ next }
    if(is.na(graph[r, "Combined_FDR"])){ next }
    
    if(graph[r, dataset] == 1){
      if(graph[r, "Combined_FDR"] == 1){
        both <- both + 1
      } else {
        in_dataset <- in_dataset + 1
      }
    } else {
      if(graph[r, "Combined_FDR"] == 1){
        in_meta <- in_meta + 1
      }
    }
  }
  graph_real <- rbind(graph_real, data.frame(dataset=dataset, category="metaanalysis", value=both))
  graph_real <- rbind(graph_real, data.frame(dataset=dataset, category="dataset", value=in_dataset))
  graph_real <- rbind(graph_real, data.frame(dataset=dataset, category="insignificant", value=in_meta))
}

graph_real <- rbind(graph_real,data.frame(dataset="metaanalysis", category="metaanalysis", value=sum(na.omit(as.numeric(graph$Combined_FDR)))))

graph_real$dataset <- gsub("FDR_", "", graph_real$dataset)


pdf("./figures/figure4.pdf", width = 10, height = 10)
ggplot(graph_real, aes(dataset, value)) + 
  geom_bar(aes(fill = reorder(graph_real$category, rev(graph_real$value))), stat = "identity", width = 0.8) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20), axis.text.y = element_text(size = 20), axis.title=element_text(size=30)) +
  guides(fill=guide_legend(title="", keywidth=0.1, keyheight=0.5, default.unit="inch")) +
  scale_x_discrete(limits = c("gCSI", "CCLE", "GDSC1000", "CTRPv2", "metaanalysis")) + 
  #scale_y_log10() + 
  labs(x = NULL, y = "number of drug tissue interactions") +
  theme(legend.text=element_text(size=30))
dev.off()

### Circos Files Generator

edges <- rbind(sig_adj_false[, c("Tissue", "Drug")], sig_adj_true[, c("Tissue", "Drug")])
edges$value <- 2
edges <- unique(edges)
rownames(edges) <- paste(edges$Tissue, edges$Drug, sep = "_")

clinical_indications <- readRDS("./clinical indications.RDS")

for(app in 1:nrow(clinical_indications)){
  if(paste(clinical_indications$variable[app], clinical_indications$drug[app], sep = "_") %in% rownames(edges)){
    edges[paste(clinical_indications$variable[app], clinical_indications$drug[app], sep = "_"), "value"] <- 3
  }else {
    edges <- rbind(edges, data.frame(Tissue=clinical_indications$variable[app], Drug=clinical_indications$drug[app], value=1))
  }
}

colnames(edges) <- c("tissue", "drug", "category")
edges <- edges[edges$drug %in% intersect(unique(c(sig_adj_true$Drug, sig_adj_false$Drug)), clinical_indications$drug),]

f <- file("karyotype.txt")
output <- c()
for(t in unique(edges$tissue))
{
  output <- c(output, paste("chr","-", t, t, 0, 100, "ylorrd-6-seq-2", sep = "\t"))
}
for(d in unique(edges$drug))
{
  output <- c(output, paste("chr", "-", d, d, 0, 100, "blues-6-seq-3", sep = "\t"))
}
writeLines(output, f, useBytes = TRUE)
close(f)


output <- c()
f <- file("edges.txt")
for(e in 1:nrow(edges))
{
  if(edges[e, "category"] == 1)
  {
    n <- paste(edges[e, "tissue"], 1, 25, edges[e, "drug"], 1, 25, "color=green_a3", sep = "\t")
  }
  else if(edges[e, "category"] == 2)
  {
    n <- paste(edges[e, "tissue"], 26, 50, edges[e, "drug"], 26, 50, "color=red_a3", sep = "\t")
  }
  else if(edges[e, "category"] == 3)
  {
    n <- paste(edges[e, "tissue"], 51, 75, edges[e, "drug"], 51, 75, "color=vdblue,z=80", sep = "\t")
  }
  else if(edges[e, "category"] == 4)
  {
    n <- paste(edges[e, "tissue"], 76, 100, edges[e, "drug"], 76, 100, "color=vlred_a3", sep = "\t")
  }
  else if(edges[e, "category"] == 5)
  {
    n <- paste(edges[e, "tissue"], 26, 50, edges[e, "drug"], 26, 50, "color=greys-5-seq-3", sep = "\t")
  }
  output <- c(n, output)
}

writeLines(output, f)
close(f)

### Jaeger Overlap (Supplementary Figure)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

download.file(url = "https://static-content.springer.com/esm/art%3A10.1186%2Fs12943-015-0312-6/MediaObjects/12943_2015_312_MOESM2_ESM.xlsx",
              destfile = "./Output/12943_2015_312_MOESM2_ESM.xlsx")
breast <- read.xlsx("./Output/12943_2015_312_MOESM2_ESM.xlsx", sheetName = "Breast cancer")
colorectal <- read.xlsx("./Output/12943_2015_312_MOESM2_ESM.xlsx", sheetName = "Colorectal cancer")
prostate <- read.xlsx("./Output/12943_2015_312_MOESM2_ESM.xlsx", sheetName = "Prostate cancer")

breast <- breast[, "Drug.name"]
colorectal <- colorectal[, "Drug.name"]
prostate <- prostate[, "Drug.name"]

meta_false <- p_adj_false[!is.na(p_adj_false$Combined_FDR), ]
meta_true <- p_adj_true[!is.na(p_adj_true$Combined_FDR), ]

clinical_indications <- readRDS("./clinical indications.RDS")

list_of_drugs <- unique(meta_false$Drug, meta_true$Drug)

#####breast
jaeger <- intersect(toupper(gsub(badchars, "", list_of_drugs)), toupper(gsub(badchars, "", breast)))
wrdm <- clinical_indications[clinical_indications$variable == "breast", "drug"]
wrdm <- toupper(gsub(badchars, "", wrdm))

vitro_p1 <- sig_adj_false[sig_adj_false$Tissue == "breast", "Drug"]
vitro_p2 <- sig_adj_true[sig_adj_true$Tissue == "breast", "Drug"]

vitro <- toupper(gsub(badchars, "", unique(c(vitro_p1, vitro_p2))))


pdf("./figures/breast.pdf", width = 20, height = 20)
draw.triple.venn(
  area1 = length(jaeger),
  area2 = length(vitro),
  area3 = length(wrdm),
  n12 = length(intersect(jaeger, vitro)),
  n13 = length(intersect(jaeger, wrdm)),
  n23 = length(intersect(vitro, wrdm)),
  n123 = length(intersect(intersect(vitro, wrdm), jaeger)),
  category = c("jaeger et al", "in vitro", "wordmined"),
  cex = 2,
  cat.cex = 2
)
dev.off()


###prostate
jaeger <- intersect(toupper(gsub(badchars, "", list_of_drugs)), toupper(gsub(badchars, "", prostate)))
wrdm <- clinical_indications[clinical_indications$variable == "prostate", "drug"]
wrdm <- toupper(gsub(badchars, "", wrdm))



vitro_p1 <- sig_adj_false[sig_adj_false$Tissue == "prostate", "Drug"]
vitro_p2 <- sig_adj_true[sig_adj_true$Tissue == "prostate", "Drug"]

vitro <- toupper(gsub(badchars, "", unique(c(vitro_p1, vitro_p2))))


pdf("./figures/prostate.pdf", width = 20, height = 20)
draw.triple.venn(
  area1 = length(jaeger),
  area2 = length(vitro),
  area3 = length(wrdm),
  n12 = length(intersect(jaeger, vitro)),
  n13 = length(intersect(jaeger, wrdm)),
  n23 = length(intersect(vitro, wrdm)),
  n123 = length(intersect(intersect(vitro, wrdm), jaeger)),
  category = c("jaeger et al", "in vitro", "wordmined"),
  cex = 2,
  cat.cex = 2
)
dev.off()



###colorectal

jaeger <- intersect(toupper(gsub(badchars, "", list_of_drugs)), toupper(gsub(badchars, "", colorectal)))
wrdm <- clinical_indications[clinical_indications$variable == "large_intestine", "drug"]
wrdm <- toupper(gsub(badchars, "", wrdm))

vitro_p1 <- sig_adj_false[sig_adj_false$Tissue == "large_intestine", "Drug"]
vitro_p2 <- sig_adj_true[sig_adj_true$Tissue == "large_intestine", "Drug"]

vitro <- toupper(gsub(badchars, "", unique(c(vitro_p1, vitro_p2))))


pdf("./figures/colorectal.pdf", width = 20, height = 20)
draw.triple.venn(
  area1 = length(jaeger),
  area2 = length(vitro),
  area3 = length(wrdm),
  n12 = length(intersect(jaeger, vitro)),
  n13 = length(intersect(jaeger, wrdm)),
  n23 = length(intersect(vitro, wrdm)),
  n123 = length(intersect(intersect(vitro, wrdm), jaeger)),
  category = c("jaeger et al", "in vitro", "wordmined"),
  cex = 2,
  cat.cex = 2
)
dev.off()

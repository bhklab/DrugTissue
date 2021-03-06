.libPaths("/rlibs/")

if(!file.exists("./temp/wordmine.RData"))
{
  source("./wordmine.R")
} else {
  load("./temp/wordmine.RData")
}

if(!file.exists("./temp/combined1.RData"))
{
  source("./gsea_with_AUC.R")
} else {
  load("./temp/combined1.RData")
}

library(VennDiagram)

path.input <- "~/Downloads"
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
breast <- xlsx::read.xlsx(paste(path.input, "12943_2015_312_MOESM2_ESM.xlsx", sep = "/"), sheetName = "Breast cancer")
colorectal <- xlsx::read.xlsx(paste(path.input, "12943_2015_312_MOESM2_ESM.xlsx", sep = "/"), sheetName = "Colorectal cancer")
prostate <- xlsx::read.xlsx(paste(path.input, "12943_2015_312_MOESM2_ESM.xlsx", sep = "/"), sheetName = "Prostate cancer")

breast <- breast[, "Drug.name"]
colorectal <- colorectal[, "Drug.name"]
prostate <- prostate[, "Drug.name"]

#breast venn diagram
options(stringsAsFactors = FALSE)

jaeger <- intersect(toupper(gsub(badchars, "", colnames(combined1))), toupper(gsub(badchars, "", breast)))

wrdm <- as.character()

for(a in 1:nrow(wordmine))
{
  if(grepl("breast", wordmine[a,"Var2"]) && wordmine[a, "value"] >= 100)
  {
    wrdm <- c(wrdm, wordmine[a, "Var1"])
  }
}
wrdm <- unique(wrdm)

c2 <- melt(as.matrix(combined1))
c2 <- na.omit(c2)
vitro <- as.character()
for(a in 1:nrow(c2))
{
  if(grepl("breast", c2[a,"Var1"]) && c2[a, "value"] <= 0.05)
  {
    vitro <- c(vitro, as.matrix(c2[a, "Var2"])[1,1])
  }
}
vitro <- unique(vitro)
vitro <- toupper(gsub(badchars, "", vitro))

pdf("breast.pdf", width = 20, height = 20)
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

#prostate
jaeger <- intersect(toupper(gsub(badchars, "", colnames(combined1))), toupper(gsub(badchars, "", prostate)))

wrdm <- as.character()

for(a in 1:nrow(wordmine))
{
  if(grepl("prostate", wordmine[a,"Var2"]) && wordmine[a, "value"] >= 100)
  {
    wrdm <- c(wrdm, wordmine[a, "Var1"])
  }
}
wrdm <- unique(wrdm)

c2 <- melt(as.matrix(combined1))
c2 <- na.omit(c2)
vitro <- as.character()
for(a in 1:nrow(c2))
{
  if(grepl("prostate", c2[a,"Var1"]) && c2[a, "value"] <= 0.05)
  {
    vitro <- c(vitro, as.matrix(c2[a, "Var2"])[1,1])
  }
}
vitro <- unique(vitro)
vitro <- toupper(gsub(badchars, "", vitro))

pdf("prostate.pdf", width = 20, height = 20)
draw.triple.venn(
  area1 = length(jaeger),
  area2 = length(vitro),
  area3 = length(wrdm),
  n12 = length(intersect(jaeger, vitro)),
  n13 = length(intersect(jaeger, wrdm)),
  n23 = length(intersect(vitro, wrdm)),
  n123 = length(intersect(intersect(vitro, wrdm), jaeger)),
  category = c("jaeger et al", "in vitro", "wordmined"),
  cat.cex = 2,
  cex = 2
)
dev.off()

#colorectal
jaeger <- intersect(toupper(gsub(badchars, "", colnames(combined1))), toupper(gsub(badchars, "", colorectal)))

wrdm <- as.character()

for(a in 1:nrow(wordmine))
{
  if(grepl("colorectal", wordmine[a,"Var2"]) && wordmine[a, "value"] >= 100)
  {
    wrdm <- c(wrdm, as.matrix(wordmine[a, "Var1"])[1,1])
  }
  if(grepl("gastrointestinal", wordmine[a,"Var2"]) && wordmine[a, "value"] >= 100)
  {
    wrdm <- c(wrdm, as.matrix(wordmine[a, "Var1"])[1,1])
  }
}
wrdm <- unique(wrdm)

c2 <- melt(as.matrix(combined1))
c2 <- na.omit(c2)
vitro <- as.character()
for(a in 1:nrow(c2))
{
  if(grepl("large_intestine", c2[a,"Var1"]) && c2[a, "value"] <= 0.05)
  {
    vitro <- c(vitro, as.matrix(c2[a, "Var2"])[1,1])
  }
}
vitro <- unique(vitro)
vitro <- toupper(gsub(badchars, "", vitro))

pdf("colorectal.pdf", width = 20, height = 20)
draw.triple.venn(
  area1 = length(jaeger),
  area2 = length(vitro),
  area3 = length(wrdm),
  n12 = length(intersect(jaeger, vitro)),
  n13 = length(intersect(jaeger, wrdm)),
  n23 = length(intersect(vitro, wrdm)),
  n123 = length(intersect(intersect(vitro, wrdm), jaeger)),
  category = c("jaeger et al", "in vitro", "wordmined"),
  cat.cex = 2,
  cex = 2,
  cex.prop = "log"
)
dev.off()
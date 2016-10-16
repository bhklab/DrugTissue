library(VennDiagram)
GDSC1000 <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")
CCLE <- downloadPSet("CCLE")
CTRPv2 <- downloadPSet("CTRPv2") 

drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="N-{3-Chloro-4-[(3-fluorobenzyl)oxy]phenyl}-6-[5-({[2-(methylsulfonyl)ethyl]amino}methyl)-2-furyl]-4-quinazolinamine", rep="lapatinib", fixed=TRUE)

pdf("figure2drugs.pdf", width = 100, height = 10, paper = "USr")
#figure 2 drugs
draw.quad.venn(
  area1 = nrow(CCLE@drug),
  area2 = nrow(GDSC1000@drug),
  area3 = nrow(CTRPv2@drug),
  area4 = nrow(gCSI@drug),
  n12 = length(intersect(rownames(CCLE@drug), rownames(GDSC1000@drug))),
  n13 = length(intersect(rownames(CCLE@drug), rownames(CTRPv2@drug))),
  n14 = length(intersect(rownames(CCLE@drug), rownames(gCSI@drug))),
  n23 = length(intersect(rownames(GDSC1000@drug), rownames(CTRPv2@drug))),
  n24 = length(intersect(rownames(GDSC1000@drug), rownames(gCSI@drug))),
  n34 = length(intersect(rownames(CTRPv2@drug), rownames(gCSI@drug))),
  n123 = length(intersect(intersect(rownames(CCLE@drug), rownames(GDSC1000@drug)), rownames(CTRPv2@drug))),
  n124 = length(intersect(intersect(rownames(CCLE@drug), rownames(GDSC1000@drug)), rownames(gCSI@drug))),
  n134 = length(intersect(intersect(rownames(CCLE@drug), rownames(CTRPv2@drug)), rownames(gCSI@drug))),
  n234 = length(intersect(intersect(rownames(GDSC1000@drug), rownames(CTRPv2@drug)), rownames(gCSI@drug))), 
  n1234 = length(intersect(intersect(rownames(CCLE@drug), rownames(GDSC1000@drug)), intersect(rownames(CTRPv2@drug), rownames(gCSI@drug)))), 
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

pdf("figure2cells.pdf", width = 100, height = 10, paper = "USr")
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



pdf("figure2tissue.pdf", width = 100, height = 10, paper = "USr")
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
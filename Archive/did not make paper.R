library(VennDiagram)

cclecell <- CCLE@cell
gdsc1000cell <- GDSC1000@cell
ctrpv2cell <- CTRPv2@cell
gcsicell <- gCSI@cell

gdsc1000cell$tissueid <- gsub("esophagus", "oesophagus", gdsc1000cell$tissueid)


cclecell <- cclecell[cclecell$tissueid != "haematopoietic_and_lymphoid_tissue", ]
gdsc1000cell <- gdsc1000cell[gdsc1000cell$tissueid != "haematopoietic_and_lymphoid_tissue", ]
ctrpv2cell <- ctrpv2cell[ctrpv2cell$tissueid != "haematopoietic_and_lymphoid_tissue", ]
gcsicell <- gcsicell[gcsicell$tissueid != "haematopoietic_and_lymphoid_tissue", ]

pdf("supp3cells.pdf", width = 100, height = 10, paper = "USr")
#figure 2 cells
draw.quad.venn(
  area1 = nrow(cclecell),
  area2 = nrow(gdsc1000cell),
  area3 = nrow(ctrpv2cell),
  area4 = nrow(gcsicell),
  n12 = length(intersect(cclecell$cellid, gdsc1000cell$cellid)),
  n13 = length(intersect(cclecell$cellid, ctrpv2cell$cellid)),
  n14 = length(intersect(cclecell$cellid, gcsicell$unique.id)),
  n23 = length(intersect(gdsc1000cell$cellid, ctrpv2cell$cellid)),
  n24 = length(intersect(gdsc1000cell$cellid, gcsicell$unique.id)),
  n34 = length(intersect(ctrpv2cell$cellid, gcsicell$unique.id)),
  n123 = length(intersect(intersect(cclecell$cellid, gdsc1000cell$cellid), ctrpv2cell$cellid)),
  n124 = length(intersect(intersect(cclecell$cellid, gdsc1000cell$cellid), gcsicell$unique.id)),
  n134 = length(intersect(intersect(cclecell$cellid, ctrpv2cell$cellid), gcsicell$unique.id)),
  n234 = length(intersect(intersect(gdsc1000cell$cellid, ctrpv2cell$cellid), gcsicell$unique.id)), 
  n1234 = length(intersect(intersect(cclecell$cellid, gdsc1000cell$cellid), intersect(ctrpv2cell$cellid, gcsicell$unique.id))), 
  category = c("CCLE", "GDSC1000", "CTRPv2", "gCSI"),
  cat.cex = 2, 
  cex = 2,
  fill = c("#F8766D", "#B79F00", "#00BFC4", "#C77CFF"),
  fontfamily = rep("Helvetica", 15),
  cat.fontfamily = rep("Helvetica", 4)
)
dev.off()

pdf("supp3tissue.pdf", width = 100, height = 10, paper = "USr")
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
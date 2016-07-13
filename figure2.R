library(VennDiagram)

pdf("figure2drugs.pdf")
#figure 2 drugs
draw.quad.venn(
  area1 = nrow(CCLE@drug),
  area2 = nrow(GDSC@drug),
  area3 = nrow(CTRPv2@drug),
  area4 = nrow(gCSI@drug),
  n12 = length(intersect(rownames(CCLE@drug), rownames(GDSC@drug))),
  n13 = length(intersect(rownames(CCLE@drug), rownames(CTRPv2@drug))),
  n14 = length(intersect(rownames(CCLE@drug), rownames(gCSI@drug))),
  n23 = length(intersect(rownames(GDSC@drug), rownames(CTRPv2@drug))),
  n24 = length(intersect(rownames(GDSC@drug), rownames(gCSI@drug))),
  n34 = length(intersect(rownames(CTRPv2@drug), rownames(gCSI@drug))),
  n123 = length(intersect(intersect(rownames(CCLE@drug), rownames(GDSC@drug)), rownames(CTRPv2@drug))),
  n124 = length(intersect(intersect(rownames(CCLE@drug), rownames(GDSC@drug)), rownames(gCSI@drug))),
  n134 = length(intersect(intersect(rownames(CCLE@drug), rownames(CTRPv2@drug)), rownames(gCSI@drug))),
  n234 = length(intersect(intersect(rownames(GDSC@drug), rownames(CTRPv2@drug)), rownames(gCSI@drug))), 
  n1234 = length(intersect(intersect(rownames(CCLE@drug), rownames(GDSC@drug)), intersect(rownames(CTRPv2@drug), rownames(gCSI@drug)))), 
  category = c("CCLE", "GDSC", "CTRPv2", "gCSI"),
  cat.cex = 2,
  cex = 2
)
dev.off()

pdf("figure2cells.pdf")
#figure 2 cells
draw.quad.venn(
  area1 = nrow(CCLE@cell),
  area2 = nrow(GDSC@cell),
  area3 = nrow(CTRPv2@cell),
  area4 = nrow(gCSI@cell),
  n12 = length(intersect(CCLE@cell$cellid, GDSC@cell$cellid)),
  n13 = length(intersect(CCLE@cell$cellid, CTRPv2@cell$cellid)),
  n14 = length(intersect(CCLE@cell$cellid, gCSI@cell$unique.id)),
  n23 = length(intersect(GDSC@cell$cellid, CTRPv2@cell$cellid)),
  n24 = length(intersect(GDSC@cell$cellid, gCSI@cell$unique.id)),
  n34 = length(intersect(CTRPv2@cell$cellid, gCSI@cell$unique.id)),
  n123 = length(intersect(intersect(CCLE@cell$cellid, GDSC@cell$cellid), CTRPv2@cell$cellid)),
  n124 = length(intersect(intersect(CCLE@cell$cellid, GDSC@cell$cellid), gCSI@cell$unique.id)),
  n134 = length(intersect(intersect(CCLE@cell$cellid, CTRPv2@cell$cellid), gCSI@cell$unique.id)),
  n234 = length(intersect(intersect(GDSC@cell$cellid, CTRPv2@cell$cellid), gCSI@cell$unique.id)), 
  n1234 = length(intersect(intersect(CCLE@cell$cellid, GDSC@cell$cellid), intersect(CTRPv2@cell$cellid, gCSI@cell$unique.id))), 
  category = c("CCLE", "GDSC", "CTRPv2", "gCSI"),
  cat.cex = 2, 
  cex = 2
)
dev.off()

pdf("figure2tissue.pdf")
#figure 2 tissue types
draw.quad.venn(
  area1 = length(unique(na.omit(CCLE@cell$tissueid))),
  area2 = length(unique(na.omit(na.omit(GDSC@cell$tissueid)))),
  area3 = length(unique(na.omit(CTRPv2@cell$tissueid[CTRPv2@cell$tissueid != ""]))),
  area4 = length(unique(na.omit(gCSI@cell$tissueid[gCSI@cell$tissueid != ""]))),
  n12 = length(intersect(na.omit(CCLE@cell$tissueid), na.omit(GDSC@cell$tissueid))),
  n13 = length(intersect(na.omit(CCLE@cell$tissueid), CTRPv2@cell$tissueid[CTRPv2@cell$tissueid != ""])),
  n14 = length(intersect(na.omit(CCLE@cell$tissueid), gCSI@cell$tissueid[gCSI@cell$tissueid != ""])),
  n23 = length(intersect(na.omit(GDSC@cell$tissueid), CTRPv2@cell$tissueid[CTRPv2@cell$tissueid != ""])),
  n24 = length(intersect(na.omit(GDSC@cell$tissueid), gCSI@cell$tissueid[gCSI@cell$tissueid != ""])),
  n34 = length(intersect(CTRPv2@cell$tissueid[CTRPv2@cell$tissueid != ""], gCSI@cell$tissueid[gCSI@cell$tissueid != ""])),
  n123 = length(intersect(intersect(na.omit(CCLE@cell$tissueid), na.omit(GDSC@cell$tissueid)), CTRPv2@cell$tissueid[CTRPv2@cell$tissueid != ""])),
  n124 = length(intersect(intersect(na.omit(CCLE@cell$tissueid), na.omit(GDSC@cell$tissueid)), gCSI@cell$tissueid[gCSI@cell$tissueid != ""])),
  n134 = length(intersect(intersect(na.omit(CCLE@cell$tissueid), CTRPv2@cell$tissueid[CTRPv2@cell$tissueid != ""]), gCSI@cell$tissueid[gCSI@cell$tissueid != ""])),
  n234 = length(intersect(intersect(na.omit(GDSC@cell$tissueid), CTRPv2@cell$tissueid[CTRPv2@cell$tissueid != ""]), gCSI@cell$tissueid[gCSI@cell$tissueid != ""])), 
  n1234 = length(intersect(intersect(na.omit(CCLE@cell$tissueid), na.omit(GDSC@cell$tissueid)), intersect(CTRPv2@cell$tissueid[CTRPv2@cell$tissueid != ""], gCSI@cell$tissueid[gCSI@cell$tissueid != ""]))), 
  category = c("CCLE", "GDSC", "CTRPv2", "gCSI"),
  cex = 2,
  cat.cex = 2
)
dev.off()
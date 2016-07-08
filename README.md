# TissueDrug

load CTRPv2, CCLE, GDSC, and gCSI

run GSEA_with_AUC.R to get enrichment scores in one matrix

make sure to run everything in one environment of RStudio as most if not all of the files depend on the variable "combined1" made from the GSEA R file. the only exception would be the wordmine files which is made from the wordmine.R file. 

Make sure to load the XML reader from drugResultGetter.R before running wordmine.R

all figures should run independently (Ctrl - A, run) , although venn diagram generation doesnt seem to wipe the previous figure and just superimposes it

lines need to be added to write to pdf automatically 

to generate the circos plot export the edges matrix generated at the end of "cytoscape edges.R" and import into cytoscape 

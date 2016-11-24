.libPaths("/rlibs/")

library(ggplot2)
library(mgcv)
library(reshape2)
library(gridExtra)
library(grid)
library(gplots)
library(RColorBrewer)
options(stringsAsFactors=FALSE)
library(openxlsx)
library(piano)
library(survcomp)
library(PharmacoGx)

options(digits = 9)

path.data=file.path("/DrugTissue/")

cellannotation <- read.csv(file.path(path.data, "cell_annotation_all_new.csv"), sep=",", comment.char="#")
breastannotation <- read.csv(file.path(path.data, "brca_cell_lines_all.csv"), sep=",", comment.char="#")

GDSC1000 <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")
CCLE <- downloadPSet("CCLE")
CTRPv2 <- downloadPSet("CTRPv2") 

dataSets <- c(gCSI, CCLE, GDSC1000)
drugs <- list()

mw <- list()
counter <- 1
ml <- list()

#removes lymphoid tissue from dataset analysis due to oversensitivity 
paper <- TRUE


drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="N-{3-Chloro-4-[(3-fluorobenzyl)oxy]phenyl}-6-[5-({[2-(methylsulfonyl)ethyl]amino}methyl)-2-furyl]-4-quinazolinamine", rep="lapatinib", fixed=TRUE)
drugNames(CCLE) <- gsub(drugNames(CCLE), pat = "AZD6244", rep = "selumetinib", fixed = TRUE)

for(d in dataSets)
{
  
  gseav <- data.frame(matrix(ncol = length(rownames(d@drug)), nrow = length(na.omit(unique(d@cell$tissueid)))))
  rownames(gseav) <- na.omit(unique(d@cell$tissueid))
  colnames(gseav) <- rownames(d@drug)
  wt <- gseav
  
  for(dr in sort(rownames(d@drug)))
  {
    if(dr != "Bleomycin (50 uM)")  
    {
      c2 <- 1
      message(paste(d@annotation$name, dr, sep = " "))
      #generate ic50 values 
      dataTable <- data.frame(matrix(ncol = 3, nrow = length(rownames(d@cell)) + length(grep("breast", d@cell$tissueid))))
      colnames(dataTable) <- c("cell_line", "unique.tissueid", dr)
      for(x in rownames(d@cell))
      {
        if(x %in% rownames(d@sensitivity$n) && dr %in% colnames(d@sensitivity$n))
        {
          if(d@sensitivity$n[x, dr] >= 1)
          {
            dataTable[c2, "cell_line"] <- x
            if(d@annotation$name == "GDSC1000" && dr == "Bleomycin")
            {
              dataTable[c2, dr] <- median(c(d@sensitivity$profiles[rownames(d@sensitivity$info)[intersect(which(x == d@sensitivity$info$cellid), which(dr == d@sensitivity$info$drugid))], "auc_recomputed"], 
                                            d@sensitivity$profiles[rownames(d@sensitivity$info)[intersect(which(x == d@sensitivity$info$cellid), which("Bleomycin (50 uM)" == d@sensitivity$info$drugid))], "auc_recomputed"])
                                          )
              
            }
            else
            {
              dataTable[c2, dr] <- median(d@sensitivity$profiles[rownames(d@sensitivity$info)[intersect(which(x == d@sensitivity$info$cellid), which(dr == d@sensitivity$info$drugid))], "auc_recomputed"])
            }
            #dataTable[c2, dr] <- median(d@sensitivity$profiles[rownames(d@sensitivity$info)[intersect(which(x == d@sensitivity$info$cellid), which(dr == d@sensitivity$info$drugid))], "auc_recomputed"])
            c2 <- c2 + 1
          }
        }
      }
      
      dataTable[, "unique.tissueid"] <- cellannotation[match(dataTable$cell_line, cellannotation[, "unique.cellid"]), "unique.tissueid"]
      
      ###RESPONSE LINE PLEASE REMOVE LATER###
      if(paper)
      {
        dataTable <- dataTable[grep("lymphoid", dataTable$unique.tissueid) * -1, ]
      }
     
      ### RESPONSE LINE ###
      
      
      ndT <- dataTable
      ch <- c2
      for(x in 1:(ch-1))
      {
        if(ndT[x, "unique.tissueid"] == "breast" && ndT[x, "unique.tissueid"] != "" && !is.na(ndT[x, "unique.tissueid"]))
        {
          ndT[c2, "cell_line"] <- ndT[x, "cell_line"]
          ndT[c2, "unique.tissueid"] <- paste("breast", breastannotation[match(ndT[x, "cell_line"], breastannotation[, "unique.cellid"]), "Subtype"], sep = "_")
          ndT[c2, dr] <- ndT[x, dr]
          c2 <- c2 + 1
        }
      }
      
      ndT <- na.omit(ndT)
      ndT <- ndT[grep("breast_NA", ndT$unique.tissueid) * -1, ]
      
      for(a in unique(ndT[, "unique.tissueid"]))
      {
        if(a != "")
        {
          wt[a, dr] <- length(grep(a, ndT[, "unique.tissueid"]))
        }
        
        if(length(grep(a, ndT[, "unique.tissueid"])) < 5 || a == "")
        {
          ndT <- ndT[ndT$unique.tissueid != a, ]
        }
      }
      
      ndT <- ndT[order(ndT$unique.tissueid), ]
      for(p in grep("esophagus", ndT$unique.tissueid))
      {
        ndT[p, "unique.tissueid"] <- "oesophagus"
      }
      ndTf <- ndT
      
      ndT <- ndTf
      if(nrow(ndT) > 0)
      {
        #linear readjustment to center around 0
        #ndT[,3] <- as.numeric(ndT[,3])
        #ndT[, 3] <- ndT[,3] - (max(ndT[,3]) + min(ndT[,3]))/2
        
        #pass into piano 
        
        #ndT[,3] <- 1/ndT[,3]
        ndT[,3] <- log10(ndT[,3])
        ndT <- ndT[order(ndT[,3]),]
        
        ndT[,3] <- 1:nrow(ndT)
        
        gset <- ndT[, c(1,2)]
        gset <- loadGSC(gset)
        
        ndTc <- ndT[grep("breast_", ndT$unique.tissueid) * -1, ]
        input <- ndTc[,3]
        names(input) <- ndTc[,1]
        
        output <- runGSA(input, gsc = gset, verbose = FALSE, nPerm = 10000)
        
        #write value to table 
        output <- GSAsummaryTable(output)
        
        if(dr == "AZD6244"){dr <- "selumetinib"}

        for(x in 1:nrow(output))
        {
          if(length(output[x, "p (non-dir.)"]) == 0)
          {
            gseav[output[x, "Name"], dr] <- NA
          }
          else
          {
            if(output[x, "p (non-dir.)"] < 1e-20)
            {
              gseav[output[x, "Name"], dr] <- 1e-20
            }else
            {
              gseav[output[x, "Name"], dr] <- output[x, "p (non-dir.)"]
            }
          }
        }
      }
    }
  }
  mw[[counter]] <- wt
  names(mw)[counter] <- d@annotation$name
  ml[[counter]] <- gseav
  names(ml)[counter] <- d@annotation$name
  counter <- counter + 1
}


save(ml, file = "./amlb.RData")
save(mw, file = "./amwb.RData")




#ctrp goes here

#drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="N-{3-Chloro-4-[(3-fluorobenzyl)oxy]phenyl}-6-[5-({[2-(methylsulfonyl)ethyl]amino}methyl)-2-furyl]-4-quinazolinamine", rep="lapatinib", fixed=TRUE)

counter <- length(dataSets) + 1
gseav <- data.frame(matrix(ncol = length(rownames(CTRPv2@drug)), nrow = length(na.omit(unique(CTRPv2@cell$tissueid)))))
rownames(gseav) <- na.omit(unique(CTRPv2@cell$tissueid))
colnames(gseav) <- rownames(CTRPv2@drug)
wt <- gseav

for(dr in sort(rownames(CTRPv2@drug)))
{

  message(paste("CTRPv2", dr, sep = " "))
  
  dataTable <- CTRPv2@sensitivity$info[CTRPv2@sensitivity$info$drugid == dr, ]
  dataTable[, "unique.tissueid"] <- CTRPv2@cell[match(dataTable$cellid, CTRPv2@cell$cellid), "tissueid"]
  dataTable[, "auc"] <- CTRPv2@sensitivity$profiles[rownames(dataTable), "auc_recomputed"]
  dataTable <- na.omit(dataTable)
  
  ndT <- dataTable
  ndT <- ndT[, c(-1, -4, -5, -6)]
  ch <- nrow(dataTable)
  c2 <- ch + 1
  for(x in 1:(ch-1))
  {
    if(ndT[x, "unique.tissueid"] == "breast" && ndT[x, "unique.tissueid"] != "" && !is.na(ndT[x, "unique.tissueid"]))
    {
      new <- data.frame("cellid" = ndT[x, "cellid"], 
                        "drugid" = ndT[x, "drugid"], 
                        "unique.tissueid" = paste("breast", breastannotation[match(ndT[x, "cellid"], breastannotation$unique.cellid), "Subtype"], sep = "_"), 
                        "auc" = ndT[x, "auc"])
      ndT <- rbind(ndT, new)
    }
  }
  
  ndT <- na.omit(ndT)
  
  for(a in unique(ndT[, "unique.tissueid"]))
  {
    if(a != "")
    {
      wt[a, dr] <- length(grep(a, ndT[, "unique.tissueid"]))
    }
    
    if(length(grep(a, ndT[, "unique.tissueid"])) < 5 || a == "")
    {
      ndT <- ndT[ndT$unique.tissueid != a, ]
    }
  }
  
  ndT <- ndT[order(ndT$unique.tissueid), ]
  if(paper)
  {
    ndT <- ndT[grep("lymphoid", ndT$unique.tissueid) * -1 , ]
  }
  ndTf <- ndT
  
  ndT <- ndTf
  if(nrow(ndT) > 0)
  {
    #linear readjustment to center around 0
    #ndT[,3] <- as.numeric(ndT[,3])
    #ndT[, 3] <- ndT[,3] - (max(ndT[,3]) + min(ndT[,3]))/2
    
    #pass into piano 
    ndT <- ndT[, -2]
    ndT[,3] <- log10(ndT[,3])
    ndT <- ndT[order(ndT[,3]),]
    ndT[,3] <- 1:nrow(ndT)
    
    gset <- ndT[, c(1,2)]
    gset <- loadGSC(gset)
    ndTc <- ndT
    if(length(grep("breast_", ndT$unique.tissueid)) > 0)
    {
      ndTc <- ndT[grep("breast_", ndT$unique.tissueid) * -1, ]
    }
    
    input <- ndTc[,3]
    names(input) <- ndTc[,1]
    
    output <- runGSA(input, gsc = gset, verbose = FALSE, nPerm = 10000)
    
    #write value to table 
    output <- GSAsummaryTable(output)
    
    if(dr == "bleomycin A2"){dr <- "Bleomycin"}
    #drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="bleomycin A2", rep="Bleomycin", fixed=TRUE)
    
    if(dr == "cytarabine hydrochloride"){dr <- "Cytarabine"}
    #drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="cytarabine hydrochloride", rep="Cytarabine", fixed=TRUE)
    
    if(dr == "doxorubicin"){dr <- "Doxorubicin"}
    #drugNames(CTRPv2) <- gsub(drugNames(CTRPv2), pat="doxorubicin", rep="Doxorubicin", fixed=TRUE)
    
    for(x in 1:nrow(output))
    {
      if(length(output[x, "p (non-dir.)"]) == 0)
      {
        gseav[output[x, "Name"], dr] <- NA
      }
      else
      {
        if(output[x, "p (non-dir.)"] < 1e-20)
        {
          gseav[output[x, "Name"], dr] <- 1e-20
        }else
        {
          gseav[output[x, "Name"], dr] <- output[x, "p (non-dir.)"]
        }
      }
    }
  }
}
mw[[counter]] <- wt
names(mw)[counter] <- "CTRPv2"
ml[[counter]] <- gseav
names(ml)[counter] <- "CTRPv2"
counter <- counter + 1

save(ml, file = "./amla.RData")
save(mw, file = "./amwa.RData")

druglist <- list()
celltypelist <- list()

for(j in ml)
{
  druglist <- c(druglist, colnames(j))
  celltypelist <- c(celltypelist, rownames(j))
}

mwc <- mw

combined <- data.frame(matrix(ncol = length(unique(druglist)), nrow = length(unique(celltypelist))))
colnames(combined) <- unique(druglist)
rownames(combined) <- unique(celltypelist)

#go through master table and get data from list of tables
for(i in rownames(combined))
{
  for(j in colnames(combined))
  {
    pvs <- vector()
    sw <- vector()
    for(k in 1:length(ml))
    {
      if(i %in% rownames(ml[[k]]) && j %in% colnames(ml[[k]]))
      {
        if(!is.na(ml[[k]][i,j]))
        {
          pvs <- c(pvs, ml[[k]][i,j])
          for(x in 1:length(mwc))
          {
            if(names(ml)[k] == names(mwc)[x])
            {
              #print(c(names(ml)[k], i, j, mwc[[x]][i,j]))
              if(j == "X17.AAG"){sw <- c(sw, mwc[[x]][i,"17.AAG"]) } 
              else if(j == "X681640"){sw <- c(sw, mwc[[x]][i,"681640"])}
              else{sw <- c(sw, mwc[[x]][i,j])}
            }
          }
        }
      }
    }
    
    if(length(pvs) != length(sw))
    {
      print(c(i,j, pvs, sw))
      message("==============================")
    }
    sw <- sw / sum(sw)
    combined[i, j] <- combine.test(pvs, weight = sw, method = "z.transform")
  }
}
combined1 <- combined
combined2 <- combined

combined1 <- melt(as.matrix(combined1))
combined1 <- na.omit(combined1)
combined1$value <- p.adjust(combined1$value, method = "fdr")
combined1 <- dcast(combined1, Var1 ~ Var2)
rownames(combined1) <- combined1$Var1
combined1 <- combined1[, -1]

for(a in colnames(combined2))
{
  combined2[,a] <- p.adjust(combined2[,a])
}

if(!paper)
{
  write.xlsx(combined1, file = "Supplementary_file_3.xlsx")
  save(combined1, file = "combined1.Rdata")
}else
{
  write.xlsx(combined1, file = "Supplementary_file_not_made.xlsx")
  save(combined1, file = "combined1mhl.Rdata")
}



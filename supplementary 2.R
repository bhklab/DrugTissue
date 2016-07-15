library(ggplot2)
library(mgcv)
library(reshape2)
library(grid)
library(gplots)
library(RColorBrewer)
options(stringsAsFactors=FALSE)
library(xlsx)
library(piano)
library(survcomp)

options(digits = 9)

path.data=file.path("data", "analysis")

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

cellannotation <- read.csv(file.path("~/Desktop/tissuedrug/PharmacoGx-private/inst/extdata", "cell_annotation_all.csv"), sep=",", comment.char="#")
breastannotation <- read.csv(file.path("~/Desktop/tissuedrug/PharmacoGx-private/inst/extdata", "brca_cell_lines_all.csv"), sep=",", comment.char="#")
dataSets <- c(NCI60)
drugs <- list()

drugsubset <- c(rownames(NCI60@curation$drug[NCI60@curation$drug$unique.drugid %in% colnames(combined1), ]))
drugsubset <- unique(drugsubset)


mw <- list()
counter <- 1
ml <- list()

for(d in dataSets)
{
  gseav <- data.frame(matrix(ncol = length(drugsubset), nrow = length(unique(na.omit(d@cell$tissueid)))))
  rownames(gseav) <- na.omit(unique(d@cell$tissueid))
  colnames(gseav) <- drugsubset
  wt <- gseav
  
  for(dr in drugsubset)
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
          dataTable[c2, dr] <- median(d@sensitivity$profiles[rownames(d@sensitivity$info)[intersect(which(x == d@sensitivity$info$cellid), which(dr == d@sensitivity$info$drugid))], "auc_recomputed"])
          c2 <- c2 + 1
        }
      }
    }
    
    dataTable[, "unique.tissueid"] <- cellannotation[match(dataTable$cell_line, cellannotation[, "unique.cellid"]), "unique.tissueid"]
    
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
    
    for(a in unique(ndT[, "unique.tissueid"]))
    {
      if(a != "")
      {
        wt[a, dr] <- length(grep(a, ndT[, "unique.tissueid"]))
      }
      
      if(a == "")
      {
        ndT <- ndT[ndT$unique.tissueid != a, ]
      }
      else if(length(grep(a, ndT[, "unique.tissueid"])) < 5 )
      {
        ndT <- ndT[ndT$unique.tissueid != a, ]
        gseav[a, dr] <- 1
      }
    }
    
    ndT <- ndT[order(ndT$unique.tissueid), ]
    ndTf <- ndT
    
    ndT <- ndTf
    if(nrow(ndT) > 0)
    {
      #linear readjustment to center around 0
      #ndT[,3] <- as.numeric(ndT[,3])
      #ndT[, 3] <- ndT[,3] - (max(ndT[,3]) + min(ndT[,3]))/2
      
      #pass into piano 
      
      #ndT[,3] <- 1/ndT[,3]
      ndT[,3] <- as.numeric(ndT[,3])
      ndT[,3] <- log10(ndT[,3])
      ndT <- ndT[order(ndT[,3]),]
      
      ndT[,3] <- 1:nrow(ndT)
      
      gset <- ndT[, c(1,2)]
      gset <- loadGSC(gset)
      
      if(TRUE %in% grepl("breast_", ndT$unique.tissueid))
      {
        ndTc <- ndT[grep("breast_", ndT$unique.tissueid) * -1, ]
      }
      else
      {
        ndTc <- ndT
      }
      input <- ndTc[,3]
      names(input) <- ndTc[,1]
      
      output <- runGSA(input, gsc = gset, verbose = FALSE, nPerm = 10000)
      
      #write value to table 
      output <- GSAsummaryTable(output)
      
      for(x in 1:nrow(output))
      {
        if(length(output[x, "p (non-dir.)"]) == 0)
        {
          gseav[output[x, "Name"], dr] <- NA
        }
        else
        {
          gseav[output[x, "Name"], dr] <- output[x, "p (non-dir.)"]
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

mwc <- mw

combined <- data.frame(matrix(ncol = length(drugsubset), nrow = length(unique(NCI60@cell$tissueid))))
colnames(combined) <- drugsubset
rownames(combined) <- unique(NCI60@cell$tissueid)

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
combined1n<- combined

combined1n <- melt(as.matrix(combined1n))
combined1n <- na.omit(combined1n)
combined1n$value <- p.adjust(combined1n$value, method = "fdr")
combined1n <- dcast(combined1n, Var1 ~ Var2)
rownames(combined1n) <- combined1n$Var1
combined1n <- combined1n[, -1]

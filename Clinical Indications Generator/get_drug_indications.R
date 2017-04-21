get_drug_indications <- function(list_of_drugs){
  
  library(XML)
  library(dplyr)
  library(plyr)
  library(reshape2)
  
  options(stringsAsFactors = FALSE)
  
  drugbank <- read.csv("~/Documents/PDX-analysis/uniprot links.csv")
  #cellannotation <- read.csv(file.path("~/Documents/DrugTissue/", "cell_annotation_all_new.csv"), sep=",", comment.char="#") 
  #tissues <- cellannotation[, grep("tissue", colnames(cellannotation))]
  drugbank <- unique(drugbank[, c("DrugBank.ID", "Name")])
  list_of_found_drugs <- c()
  output <- data.frame()
  for(test_drug in list_of_drugs){
    
    drugbank_id <- drugbank[which(toupper(test_drug) == toupper(drugbank$Name)), "DrugBank.ID"]
    
    if(length(drugbank_id) != 0){
      
      list_of_found_drugs <- c(list_of_found_drugs, test_drug)
      
      message(test_drug)
      link_to_parse <- paste("https://www.drugbank.ca/drugs/", drugbank_id, sep = "")
      drug_data <- readLines(link_to_parse, warn = FALSE)
      parsed_data <- htmlParse(drug_data, asTree = TRUE)
      
      raw <- xpathSApply(parsed_data, "//li/a[contains(@href, 'indications')]", xmlValue)
      raw_p2 <-  xpathSApply(parsed_data, "//tr[./th = 'Indication']/td", xmlValue)
      
      drug_identifiers <- data.frame(matrix(ncol = length(raw) + 1, nrow = 1))
      colnames(drug_identifiers) <- c("drug", as.vector(raw))
      drug_identifiers[1,] <- TRUE
      drug_identifiers[1,1] <- test_drug
      rownames(drug_identifiers)[1] <- test_drug
      output <- rbind.fill(output, drug_identifiers)
    }
  }
  return(output)
}


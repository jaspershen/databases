masstools::setwd_project()
rm(list= ls())
setwd("raw_data/KEGG/")

library(metid)
library(tidyverse)

library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)

rm(list = ls())

###to get KEGG database
library(KEGGgraph)
library(KEGGREST)
library(KEGGlincs)
library(tidyverse)

compound_ID <-
  keggList(database = "compound") %>%
  names() %>%
  unique() %>%
  stringr::str_replace_all(., "cpd:", "")

# kegg_compound_database <-
#   pbapply::pblapply(compound_ID, function(x){
#     KEGGREST::keggGet(dbentries = x)[[1]]
#   })
# 
# save(kegg_compound_database, file = "kegg_compound_database")

load("kegg_compound_database")

kegg_metabolite =
  kegg_compound_database %>%
  purrr::map(function(x) {
    cat(x$ENTRY, " ")
    KEGG.ID = x$ENTRY
    x$NAME <- stringr::str_replace(x$NAME, "\\;$", "")
    Compound.name = paste(x$NAME, collapse = "{}")
    Formula = x$FORMULA
    if (is.null(x$FORMULA)) {
      Formula = NA
    }
    mz = as.numeric(x$EXACT_MASS)
    if (is.null(x$EXACT_MASS)) {
      mz = NA
    }
    CAS.ID = stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
      stringr::str_trim(side = "both")
    
    PubChem.ID = stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
      stringr::str_trim(side = "both")
    
    if (length(CAS.ID) == 0) {
      CAS.ID = NA
    }
    
    if (length(PubChem.ID) == 0) {
      PubChem.ID = NA
    }
    
    From_human = TRUE
    REMARK <- x$REMARK
    if(is.null(REMARK)){
      From_drug <- FALSE
      KEGG_drug.ID <- NA
    }else{
      KEGG_drug.ID <- 
        paste(stringr::str_extract_all(REMARK, "D[0-9]{5,6}")[[1]], collapse = "{}")
      
      if(length(KEGG_drug.ID) == 0){
        From_drug <- FALSE
        KEGG_drug.ID <- NA
      }else{
        From_drug <- TRUE
        KEGG_drug.ID <- KEGG_drug.ID
      }   
    }
    
    data.frame(
      Lab.ID = KEGG.ID,
      Compound.name,
      Formula,
      mz,
      CAS.ID,
      HMDB.ID = NA,
      KEGG.ID,
      PubChem.ID,
      From_human = From_human,
      From_drug = From_drug,
      KEGG_drug.ID = KEGG_drug.ID 
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

kegg_metabolite <-
  kegg_metabolite %>%
  dplyr::filter(!is.na(mz)) %>%
  dplyr::mutate(synonym = Compound.name) %>%
  dplyr::mutate(
    RT = NA,
    mz.pos = NA,
    mz.neg = NA,
    Submitter = "KEGG"
  ) %>%
  dplyr::select(
    Lab.ID,
    Compound.name,
    mz,
    RT,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Formula,
    mz.pos,
    mz.neg,
    Submitter,
    everything()
  )

kegg_metabolite$Compound.name =
  kegg_metabolite$Compound.name %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist() %>%
  stringr::str_replace(pattern = ";", "")

kegg_metabolite$synonym


save(kegg_metabolite, file = "kegg_metabolite")







drug_ID <-
  keggList(database = "drug") %>%
  names() %>%
  unique() %>%
  stringr::str_replace_all(., "dr:", "")

# kegg_drug_database <-
#   pbapply::pblapply(drug_ID, function(x){
#     KEGGREST::keggGet(dbentries = x)[[1]]
#   })
# 
# save(kegg_drug_database, file = "kegg_drug_database")

load("kegg_compound_database")

kegg_metabolite =
  kegg_compound_database %>%
  purrr::map(function(x) {
    cat(x$ENTRY, " ")
    KEGG.ID = x$ENTRY
    x$NAME <- stringr::str_replace(x$NAME, "\\;$", "")
    Compound.name = paste(x$NAME, collapse = "{}")
    Formula = x$FORMULA
    if (is.null(x$FORMULA)) {
      Formula = NA
    }
    mz = as.numeric(x$EXACT_MASS)
    if (is.null(x$EXACT_MASS)) {
      mz = NA
    }
    CAS.ID = stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
      stringr::str_trim(side = "both")
    
    PubChem.ID = stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
      stringr::str_trim(side = "both")
    
    if (length(CAS.ID) == 0) {
      CAS.ID = NA
    }
    
    if (length(PubChem.ID) == 0) {
      PubChem.ID = NA
    }
    
    From_human = TRUE
    REMARK <- x$REMARK
    if(is.null(REMARK)){
      From_drug <- FALSE
      KEGG_drug.ID <- NA
    }else{
      KEGG_drug.ID <- 
        paste(stringr::str_extract_all(REMARK, "D[0-9]{5,6}")[[1]], collapse = "{}")
      
      if(length(KEGG_drug.ID) == 0){
        From_drug <- FALSE
        KEGG_drug.ID <- NA
      }else{
        From_drug <- TRUE
        KEGG_drug.ID <- KEGG_drug.ID
      }   
    }
    
    data.frame(
      Lab.ID = KEGG.ID,
      Compound.name,
      Formula,
      mz,
      CAS.ID,
      HMDB.ID = NA,
      KEGG.ID,
      PubChem.ID,
      From_human = From_human,
      From_drug = From_drug,
      KEGG_drug.ID = KEGG_drug.ID 
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

kegg_metabolite <-
  kegg_metabolite %>%
  dplyr::filter(!is.na(mz)) %>%
  dplyr::mutate(synonym = Compound.name) %>%
  dplyr::mutate(
    RT = NA,
    mz.pos = NA,
    mz.neg = NA,
    Submitter = "KEGG"
  ) %>%
  dplyr::select(
    Lab.ID,
    Compound.name,
    mz,
    RT,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Formula,
    mz.pos,
    mz.neg,
    Submitter,
    everything()
  )

kegg_metabolite$Compound.name =
  kegg_metabolite$Compound.name %>%
  stringr::str_split(pattern = "\\{\\}") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist() %>%
  stringr::str_replace(pattern = ";", "")

kegg_metabolite$synonym


save(kegg_metabolite, file = "kegg_metabolite")


# keggMS1Database_1.0 =
#   construct_database(
#     path = ".",
#     version = "1.0",
#     metabolite.info.name = "kegg.xlsx",
#     source = "KEGG",
#     link = "https://www.genome.jp/kegg/compound/",
#     creater = "Xiaotao Shen",
#     email = "shenxt@stanford.edu",
#     rt = FALSE,
#     threads = 3
#   )
# 
# save(keggMS1Database_1.0, file = "keggMS1Database_1.0")

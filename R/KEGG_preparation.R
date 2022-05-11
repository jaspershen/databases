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
library(KEGGREST)
library(tidyverse)

# # compound_ID <-
# #   keggList(database = "compound") %>%
# #   names() %>%
# #   unique() %>%
# #   stringr::str_replace_all(., "cpd:", "")
# #
# # kegg_compound_database <-
# #   pbapply::pblapply(compound_ID, function(x){
# #     KEGGREST::keggGet(dbentries = x)[[1]]
# #   })
# #
# # save(kegg_compound_database, file = "kegg_compound_database")
# 
# load("kegg_compound_database")
# 
# kegg_compound_database %>%
#   lapply(function(x){
#     stringr::str_split(x$DBLINK, ": ") %>%
#       lapply(function(x){
#         x[1]
#       }) %>%
#       unlist()
#   }) %>%
#   unlist() %>%
#   unique()
# 
# kegg_metabolite =
#   kegg_compound_database %>%
#   purrr::map(function(x) {
#     cat(x$ENTRY, " ")
#     KEGG.ID = x$ENTRY
#     x$NAME <- stringr::str_replace(x$NAME, "\\;$", "")
#     Compound.name = paste(x$NAME, collapse = "{}")
#     Formula = x$FORMULA
#     if (is.null(x$FORMULA)) {
#       Formula = NA
#     }
#     mz = as.numeric(x$EXACT_MASS)
#     if (is.null(x$EXACT_MASS)) {
#       mz = NA
#     }
# 
#     CAS.ID = stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     PUBCHEM.ID = stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     CHEBI.ID <-
#       stringr::str_replace(grep("ChEBI", x$DBLINKS, value = TRUE), "ChEBI: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     CHEMBL.ID <-
#       stringr::str_replace(grep("ChEMBL", x$DBLINKS, value = TRUE), "ChEMBL: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     LIPIDMAPS.ID <-
#       stringr::str_replace(grep("LIPIDMAPS", x$DBLINKS, value = TRUE), "LIPIDMAPS: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     LIPIDBANK.ID <-
#       stringr::str_replace(grep("LipidBank", x$DBLINKS, value = TRUE), "LipidBank: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     DRUGBANK.ID <-
#       stringr::str_replace(grep("DrugBank", x$DBLINKS, value = TRUE), "DrugBank: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     if (length(CAS.ID) == 0) {
#       CAS.ID = NA
#     }
# 
#     if (length(PUBCHEM.ID) == 0) {
#       PUBCHEM.ID = NA
#     }
# 
#     if (length(CHEBI.ID) == 0) {
#       CHEBI.ID = NA
#     }
# 
#     if (length(CHEMBL.ID) == 0) {
#       CHEMBL.ID = NA
#     }
# 
#     if (length(LIPIDMAPS.ID) == 0) {
#       LIPIDMAPS.ID = NA
#     }
# 
#     if (length(LIPIDBANK.ID) == 0) {
#       LIPIDBANK.ID = NA
#     }
# 
#     if (length(DRUGBANK.ID) == 0) {
#       DRUGBANK.ID = NA
#     }
# 
#     From_human = "Yes"
#     REMARK <- x$REMARK
#     if(is.null(REMARK)){
#       From_drug <- "No"
#       KEGG_DRUG.ID <- NA
#     }else{
#       KEGG_DRUG.ID <-
#         paste(stringr::str_extract_all(REMARK, "D[0-9]{5,6}")[[1]], collapse = "{}")
# 
#       if(length(KEGG_DRUG.ID) == 0){
#         From_drug <- "No"
#         KEGG_DRUG.ID <- NA
#       }else{
#         From_drug <- "Yes"
#         KEGG_DRUG.ID <- KEGG_DRUG.ID
#       }
#     }
# 
#     data.frame(
#       Lab.ID = KEGG.ID,
#       Compound.name,
#       Formula,
#       mz,
#       CAS.ID,
#       HMDB.ID = NA,
#       KEGG.ID,
#       PUBCHEM.ID,
#       CHEBI.ID,
#       CHEMBL.ID,
#       LIPIDMAPS.ID,
#       LIPIDBANK.ID,
#       DRUGBANK.ID = DRUGBANK.ID,
#       From_human = From_human,
#       From_drug = From_drug,
#       KEGG_DRUG.ID = KEGG_DRUG.ID
#     )
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# kegg_metabolite <-
#   kegg_metabolite %>%
#   dplyr::filter(!is.na(mz)) %>%
#   dplyr::mutate(Synonyms = Compound.name) %>%
#   dplyr::mutate(
#     RT = NA,
#     mz.pos = NA,
#     mz.neg = NA,
#     Submitter = "KEGG"
#   ) %>%
#   dplyr::select(
#     Lab.ID,
#     Compound.name,
#     mz,
#     RT,
#     CAS.ID,
#     HMDB.ID,
#     KEGG.ID,
#     Formula,
#     mz.pos,
#     mz.neg,
#     Submitter,
#     everything()
#   )
# 
# kegg_metabolite$Compound.name =
#   kegg_metabolite$Compound.name %>%
#   stringr::str_split(pattern = "\\{\\}") %>%
#   purrr::map(function(x) {
#     x[1]
#   }) %>%
#   unlist() %>%
#   stringr::str_replace(pattern = ";", "")
# 
# 
# save(kegg_metabolite, file = "kegg_metabolite")
load("kegg_metabolite")

# drug_ID <-
#   keggList(database = "drug") %>%
#   names() %>%
#   unique() %>%
#   stringr::str_replace_all(., "cpd:", "")
# 
# kegg_drug_database <-
#   pbapply::pblapply(drug_ID, function(x){
#     KEGGREST::keggGet(dbentries = x)[[1]]
#   })
#
# save(kegg_drug_database, file = "kegg_drug_database")

load("kegg_drug_database")

# kegg_drug_database %>%
#   lapply(function(x){
#     stringr::str_split(x$DBLINK, ": ") %>%
#       lapply(function(x){
#         x[1]
#       }) %>%
#       unlist()
#   }) %>%
#   unlist() %>%
#   unique()
# 
# kegg_drug =
#   kegg_drug_database %>%
#   purrr::map(function(x) {
#     cat(x$ENTRY, " ")
#     KEGG_DRUG.ID = x$ENTRY
#     x$NAME <- stringr::str_replace(x$NAME, "\\;$", "")
#     Compound.name = paste(x$NAME, collapse = "{}")
#     Formula = x$FORMULA
#     if (is.null(x$FORMULA)) {
#       Formula = NA
#     }
#     mz = as.numeric(x$EXACT_MASS)
#     if (is.null(x$EXACT_MASS)) {
#       mz = NA
#     }
# 
#     CAS.ID = stringr::str_replace(grep("CAS", x$DBLINKS, value = TRUE), "CAS: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     PUBCHEM.ID = stringr::str_replace(grep("PubChem", x$DBLINKS, value = TRUE), "PubChem: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     CHEBI.ID <-
#       stringr::str_replace(grep("ChEBI", x$DBLINKS, value = TRUE), "ChEBI: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     CHEMBL.ID <-
#       stringr::str_replace(grep("ChEMBL", x$DBLINKS, value = TRUE), "ChEMBL: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     LIPIDMAPS.ID <-
#       stringr::str_replace(grep("LIPIDMAPS", x$DBLINKS, value = TRUE), "LIPIDMAPS: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     LIPIDBANK.ID <-
#       stringr::str_replace(grep("LipidBank", x$DBLINKS, value = TRUE), "LipidBank: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     DRUGBANK.ID <-
#       stringr::str_replace(grep("DrugBank", x$DBLINKS, value = TRUE), "DrugBank: ", "") %>%
#       stringr::str_trim(side = "both")
# 
#     if (length(CAS.ID) == 0) {
#       CAS.ID = NA
#     }
# 
#     if (length(PUBCHEM.ID) == 0) {
#       PUBCHEM.ID = NA
#     }
# 
#     if (length(CHEBI.ID) == 0) {
#       CHEBI.ID = NA
#     }
# 
#     if (length(CHEMBL.ID) == 0) {
#       CHEMBL.ID = NA
#     }
# 
#     if (length(LIPIDMAPS.ID) == 0) {
#       LIPIDMAPS.ID = NA
#     }
# 
#     if (length(LIPIDBANK.ID) == 0) {
#       LIPIDBANK.ID = NA
#     }
# 
#     if (length(DRUGBANK.ID) == 0) {
#       DRUGBANK.ID = NA
#     }
# 
#     From_drug = "Yes"
#     REMARK <- x$REMARK
#     if(is.null(REMARK)){
#       From_human <- "No"
#       KEGG.ID <- NA
#     }else{
#       KEGG.ID <-
#         stringr::str_extract_all(REMARK, "Same as: C[0-9]{5,6}")[[1]] %>%
#         stringr::str_replace_all("Same as: ", "") %>%
#         paste(collapse = "{}")
# 
#       if(KEGG.ID == ""){
#         From_human <- "No"
#         KEGG.ID <- NA
#       }else{
#         From_human <- "Yes"
#         KEGG.ID <- KEGG.ID
#       }
#     }
# 
#     data.frame(
#       Lab.ID = KEGG_DRUG.ID,
#       Compound.name,
#       Formula,
#       mz,
#       CAS.ID,
#       HMDB.ID = NA,
#       KEGG.ID,
#       PUBCHEM.ID,
#       CHEBI.ID,
#       CHEMBL.ID,
#       LIPIDMAPS.ID,
#       LIPIDBANK.ID,
#       DRUGBANK.ID = DRUGBANK.ID,
#       From_human = From_human,
#       From_drug = From_drug,
#       KEGG_DRUG.ID = KEGG_DRUG.ID
#     )
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# kegg_drug <-
#   kegg_drug %>%
#   dplyr::filter(!is.na(mz)) %>%
#   dplyr::mutate(Synonyms = Compound.name) %>%
#   dplyr::mutate(
#     RT = NA,
#     mz.pos = NA,
#     mz.neg = NA,
#     Submitter = "KEGG"
#   ) %>%
#   dplyr::select(
#     Lab.ID,
#     Compound.name,
#     mz,
#     RT,
#     CAS.ID,
#     HMDB.ID,
#     KEGG.ID,
#     Formula,
#     mz.pos,
#     mz.neg,
#     Submitter,
#     everything()
#   )
# 
# kegg_drug$Compound.name =
#   kegg_drug$Compound.name %>%
#   stringr::str_split(pattern = "\\{\\}") %>%
#   purrr::map(function(x) {
#     x[1]
#   }) %>%
#   unlist() %>%
#   stringr::str_replace(pattern = ";", "")
# 
# kegg_drug[which(kegg_drug == "", arr.ind = TRUE)] <- NA
# 
# save(kegg_drug, file = "kegg_drug")

load("kegg_drug")

dim(kegg_metabolite)
dim(kegg_drug)

kegg_metabolite[which(kegg_metabolite == "", arr.ind = TRUE)] <- NA
kegg_drug[which(kegg_drug == "", arr.ind = TRUE)] <- NA

kegg_drug.id <- 
kegg_metabolite %>% 
  dplyr::filter(!is.na(KEGG_DRUG.ID)) %>% 
  dplyr::pull(KEGG_DRUG.ID) %>% 
  stringr::str_split("\\{\\}") %>% 
  unlist() %>% 
  unique()

kegg_drug <- 
  kegg_drug %>% 
  dplyr::filter(!KEGG_DRUG.ID %in% kegg_drug.id)
  
kegg.id <- 
kegg_drug$KEGG.ID[!is.na(kegg_drug$KEGG.ID)]

kegg_metabolite <- 
kegg_metabolite %>% 
  dplyr::filter(!KEGG.ID %in% kegg.id)

kegg_metabolite$From_human
sum(kegg_metabolite$From_drug == "Yes")

sum(kegg_drug$From_human == "Yes")
sum(kegg_drug$From_drug == "Yes")


colnames(kegg_metabolite)
colnames(kegg_drug)

kegg_ms1 <- 
  rbind(kegg_metabolite,
        kegg_drug)

openxlsx::write.xlsx(kegg_ms1, file = "kegg.xlsx", asTable = TRUE)

library(metid)

kegg_ms1 =
  construct_database(
    path = ".",
    version = "2022-04-08",
    metabolite.info.name = "kegg.xlsx",
    source = "KEGG",
    link = "https://www.genome.jp/kegg",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

masstools::setwd_project()
setwd("data/KEGG")
save(kegg_ms1, file = "kegg_ms1.rda", compress = "xz")






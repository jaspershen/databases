## http://bigg.ucsd.edu/
## 
no_source()
masstools::setwd_project()
rm(list=ls())
source("R/BIGG.R")
setwd("raw_data/BIGG/")
library(jsonlite)
library(rjson)

#####download universal metabolites
metabolite_info <- 
  request_bigg_universal_metabolite_info()

# bigg_metabolite <-
#   vector(mode = "list", length = nrow(metabolite_info))
# 
# for(i in 1:length(bigg_metabolite)){
#   cat(i, " ")
#   Sys.sleep(time = "5")
#   bigg_metabolite[[i]] <-
#   request_bigg_universal_metabolite(metabolite_id = metabolite_info$bigg_id[i], return_form = "data.frame")
# }
# 
# save(bigg_metabolite, file = "bigg_metabolite")

load("bigg_metabolite")

# lapply(bigg_metabolite, function(x){
#   any(colnames(x) == "Var.8")
# }) %>% 
#   unlist() %>% 
#   which()

bigg_metabolite <- 
  bigg_metabolite %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

#####model information
# model_info <- 
#   request_bigg_model_info()
# 
# write.csv(model_info, "model_info.csv", row.names = FALSE)
# 
# model_info <- readr::read_csv("model_info_manual.csv")
# 
# model_class <- 
# bigg_metabolite$model_bigg_id %>% 
#   lapply(function(x){
#     model_class <- 
#     model_info %>% 
#       dplyr::filter(bigg_id %in% unique(stringr::str_split(x, "\\{\\}")[[1]])) %>% 
#       pull(class) %>% 
#       unique() %>% 
#       sort() %>% 
#       paste(collapse = "{}")
#     model_class
#   }) %>% 
#   unlist()
# 
# metabolite_soure <-
#   matrix("No", nrow = nrow(bigg_metabolite), ncol = 13) %>%
#   as.data.frame()
# 
# colnames(metabolite_soure) <-
#   c("From_human", "From_mammal", "From_microbiota", "From_archaea",
#     "From_bacteria", "From_fungi", "From_eukaryota", "From_food_plant",
#     "From_food", "From_plan", "From_drug", "From_environment", "From_other")
# 
# metabolite_soure
# 
# for (i in 1:length(model_class)) {
#   cat(i, " ")
#   x = model_class[i]
#   x <-
#     stringr::str_split(x, pattern = "\\{\\}")[[1]]
#   if (any(x %in% c("Archaea"))) {
#     metabolite_soure$From_archaea[i] <- "Yes"
#     metabolite_soure$From_microbiota[i] <- "Yes"
#   }
# 
#   if (any(x %in% c("Bacteria"))) {
#     metabolite_soure$From_bacteria[i] <- "Yes"
#     metabolite_soure$From_microbiota[i] <- "Yes"
#   }
# 
#   if (any(x %in% c("Fungi"))) {
#     metabolite_soure$From_fungi[i] <- "Yes"
#     metabolite_soure$From_microbiota[i] <- "Yes"
#   }
# 
#   if (any(x %in% c("Eukaryota"))) {
#     metabolite_soure$From_eukaryota[i] <- "Yes"
#     metabolite_soure$From_microbiota[i] <- "Yes"
#   }
# 
#   if (any(x %in% c("Human"))) {
#     metabolite_soure$From_human[i] <- "Yes"
#   }
# 
#   if (any(x %in% c("Mammal"))) {
#     metabolite_soure$From_mammal[i] <- "Yes"
#   }
# }
# 
# save(metabolite_soure, file = "metabolite_soure")
load("metabolite_soure")

bigg_metabolite <- 
cbind(bigg_metabolite, 
      metabolite_soure)


colnames(bigg_metabolite)

####Formula and mz
bigg_metabolite[which(bigg_metabolite == "", arr.ind = TRUE)] <- NA

bigg_metabolite %>% 
  dplyr::filter(is.na(formula)) %>% 
  pull(bigg_id)

bigg_metabolite <- 
bigg_metabolite %>% 
  dplyr::filter(!is.na(formula))


bigg_metabolite %>% 
  dplyr::filter(!is.na(KEGG)) %>% 
  dplyr::select(bigg_id, charges, formula, KEGG, HMDB)

####add formula to H or remove H based on charge

formula <- 
1:nrow(bigg_metabolite) %>% 
  purrr::map(function(idx){
    cat(idx, " ")
    charge = bigg_metabolite$charges[idx]
    formula <- bigg_metabolite$formula[idx]
    
    if(is.na(charge)){
      return(bigg_metabolite$formula[idx])
    }
    
    if(charge == "0"){
      return(bigg_metabolite$formula[idx])
    }
    
    if(length(grep("\\{\\}",charge)) == 0){
      charge <- as.numeric(charge)
      adduct = paste0(ifelse(charge > 0, "M-", "M+"), abs(charge), "H")
      formula <- masstools::sum_formula(formula = formula, adduct = adduct)
      return(formula)
    }
    
    if(length(grep("\\{\\}",charge)) > 0){
      charge <- as.numeric(stringr::str_split(charge, "\\{\\}")[[1]])
      formula <- stringr::str_split(formula, "\\{\\}")[[1]]
      if(length(charge) != length(formula)){
        return(NA)
      }
      
      if(any(charge) == 0){
        return(formula[which(charge == 0)])
      }
      
      adduct = paste0(ifelse(charge[1] > 0, "M-", "M+"), abs(charge[1]), "H")
      formula <- masstools::sum_formula(formula = formula[1], adduct = adduct)
      return(formula)
    }
  }) %>% 
  unlist()

bigg_metabolite$formula <- formula


load("../../data/KEGG/kegg_ms1.rda")

temp <- 
bigg_metabolite[,c("bigg_id","KEGG", "formula")] %>% 
  dplyr::filter(!is.na(KEGG)) %>% 
  dplyr::filter(!is.na(formula)) %>% 
  dplyr::left_join(kegg_ms1@spectra.info[,c("KEGG.ID", "Formula")],
                   by = c("KEGG" = "KEGG.ID")) %>% 
  dplyr::filter(!is.na(Formula))

temp %>% 
  dplyr::filter(Formula != formula)
  

formula[match(temp$bigg_id, bigg_metabolite$bigg_id)] <- temp$Formula

bigg_metabolite$formula <- formula


bigg_metabolite <- 
  bigg_metabolite %>% 
  dplyr::filter(!is.na(formula))

###add mz
library(Rdisop)
mz <-
  bigg_metabolite$formula %>% 
  purrr::map(function(x){
  temp <- tryCatch(Rdisop::getMass(getMolecule(x)), error = function(e)NA)
  }) %>% 
  unlist()

bigg_metabolite$mz <- mz

bigg_metabolite <- 
bigg_metabolite %>% 
  dplyr::filter(!is.na(mz))


bigg_metabolite <- 
bigg_metabolite %>% 
  dplyr::rename(BIGG.ID = bigg_id,
                Compound.name = name,
                Formula = formula,
                BIOCYC.ID = BioCyc,
                CHEBI.ID = CHEBI,
                HMDB.ID = HMDB,
                INCHIKEY.ID = InChI_Key,
                KEGG.ID = KEGG,
                METANETX.ID = MetaNetX,
                REACTOME.ID = Reactome,
                SEED.ID = SEED,
                KEGG_DRUG.ID = KEGG_Drug,
                KEGG_GLYCAN.ID = KEGG_Glycan,
                LIPIDMAPS.ID = LipidMaps
                ) %>% 
  dplyr::mutate(Lab.ID = BIGG.ID,
                CAS.ID = NA,
                RT = NA,
                mz.pos = NA,
                mz.neg = NA,
                Submitter = "BIGG") %>% 
  dplyr::select(Lab.ID,
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
                everything())

bigg_metabolite$Synonyms <- 
  bigg_metabolite$Compound.name %>% 
  stringr::str_replace_all("; ", "\\{\\}")

bigg_metabolite$Compound.name <- 
  bigg_metabolite$Synonyms %>% 
  stringr::str_split("\\{\\}") %>% 
  lapply(function(x){
    x[1]
  }) %>% 
  unlist()

bigg_metabolite$HMDB.ID


openxlsx::write.xlsx(bigg_metabolite, file = "bigg_metabolite.xlsx", asTable = TRUE)

bigg_ms1=
  construct_database(
    path = ".",
    version = "5.1.8",
    metabolite.info.name = "bigg_metabolite.xlsx",
    source = "BIGG",
    link = "http://bigg.ucsd.edu/",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

masstools::setwd_project()
setwd("data/BIGG")
save(bigg_ms1, file = "bigg_ms1.rda", compress = "xz")




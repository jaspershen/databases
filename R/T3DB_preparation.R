masstools::setwd_project()
setwd('raw_data/T3DB')
library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)
rm(list = ls())

library(rjson)

t3db_2 =
  readr::read_csv("toxins.csv")

t3db =
  jsonlite::fromJSON(txt = "toxins.json")

t3db <-
  t3db %>%
  dplyr::left_join(t3db_2[, c("T3DB ID", "Class")], by = c("title" = "T3DB ID")) %>%
  dplyr::filter(Class == "SmallMolecule")

colnames(t3db)

t3db <-
  t3db %>%
  dplyr::filter(Class == "SmallMolecule") %>%
  dplyr::select(-id) %>% 
  dplyr::rename(
    Lab.ID = title,
    PUBCHEM.ID = pubchem_id,
    Compound.name = common_name,
    HMDB.ID = hmdb_id,
    CAS.ID = cas,
    Formula = chemical_formula,
    mz = moldb_mono_mass,
    WIKIPEDIA.ID = wikipedia,
    KEGG.ID = kegg_compound_id,
    UNIPROT.ID = uniprot_id,
    OMIN.ID = omim_id,
    CHEBI.ID = chebi_id,
    BIOCYC.ID = biocyc_id,
    CTD.ID = ctd_id,
    STITCH.ID = stitch_id,
    DRUGBANK.ID = drugbank_id,
    PDB.ID = pdb_id,
    CHEMBL.ID = chembl_id,
    CHEMSPIDER.ID = chemspider_id,
    BIODB.ID = biodb_id,
    SMILES.ID = moldb_smiles,
    INCHI.ID = moldb_inchi,
    INCHIKEY.ID = moldb_inchikey, 
    Create_date = created_at,
    Updated_date = updated_at,
    Average.mass = moldb_average_mass,
    Synonyms = synonyms_list
  ) %>%
  dplyr::mutate(
    T3DB.ID = Lab.ID,
    RT = NA,
    mz.pos = NA,
    mz.neg = NA,
    Submitter = "T3DB"
  ) %>%
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

# t3db$Formula
t3db$From_environment <- "Yes"

t3db$toxicity <- 
  t3db$toxicity %>% 
  stringr::str_replace_all('\r\n', "{}")

t3db$Synonyms <-   
t3db$Synonyms %>% 
  stringr::str_replace_all('\r\n', "{}")

t3db$types <-   
  t3db$types %>% 
  lapply(function(x){
    if(nrow(x) == 0){
      return(NA)
    }
    paste(x$type_name, collapse = "{}")
  }) %>% 
  unlist()

t3db$cellular_locations <-   
  t3db$cellular_locations %>% 
  lapply(function(x){
    if(nrow(x) == 0){
      return(NA)
    }
    paste(x$name, collapse = "{}")
  }) %>% 
  unlist()

t3db <- 
t3db %>% 
  dplyr::rename(Cellular_locations = cellular_locations)


t3db$tissues <-   
  t3db$tissues %>% 
  lapply(function(x){
    if(nrow(x) == 0){
      return(NA)
    }
    paste(x$name, collapse = "{}")
  }) %>% 
  unlist()

t3db <- 
  t3db %>% 
  dplyr::rename(Tissues = tissues)

t3db$Pathway_name <-   
  t3db$pathways %>% 
  lapply(function(x){
    if(nrow(x) == 0){
      return(NA)
    }
    paste(x$name, collapse = "{}")
  }) %>% 
  unlist()

t3db$Pathway_KEGG.ID <-   
  t3db$pathways %>% 
  lapply(function(x){
    if(nrow(x) == 0){
      return(NA)
    }
    paste(x$kegg_id[!is.na(x$kegg_id)], collapse = "{}")
  }) %>% 
  unlist()

t3db$Pathway_SMPDB.ID <-   
  t3db$pathways %>% 
  lapply(function(x){
    if(nrow(x) == 0){
      return(NA)
    }
    x$kegg_id[x$kegg_id == ""] <- NA
    x$smpdb_id[x$smpdb_id == ""] <- NA
    if(all(is.na(x$smpdb_id))){
      return(NA)
    }
    paste(x$smpdb_id[!is.na(x$smpdb_id)], collapse = "{}")
  }) %>% 
  unlist()

t3db <-
  t3db %>%
  dplyr::select(-pathways)

t3db[which(t3db == "", arr.ind = TRUE)] <- NA


load("../KEGG/kegg_ms1.rda")


t3db$KEGG.ID

t3db <-
  t3db %>% dplyr::left_join(kegg_ms1@spectra.info[, c("KEGG.ID", "From_human", "From_drug")],
                            by = c("KEGG.ID"),
                            na_matches = "never")

t3db$From_drug[is.na(t3db$From_drug)] <- "No"
t3db$From_human[is.na(t3db$From_human)] <- "No"

openxlsx::write.xlsx(t3db, file = "t3db.xlsx", asTable = TRUE)

library(metid)

t3db_ms1 =
  construct_database(
    path = ".",
    version = "2022-04-08",
    metabolite.info.name = "t3db.xlsx",
    source = "T3DB",
    link = "http://www.t3db.ca/",
    creater = "Xiaotao Shen",
    email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = 3
  )

masstools::setwd_project()
setwd("data/T3DB")
save(t3db_ms1, file = "t3db_ms1.rda", compress = "xz")



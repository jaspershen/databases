##
no_source()
masstools::setwd_project()
library(dplyr)
library(ggplot2)
library(XML)
library(MetaDBparse)
rm(list = ls())

library(RJSONIO)
library(readr)

library(xml2)

source("R/foodb.R")

setwd('raw_data/foodb')

compound <-
  readr::read_csv("foodb_2020_04_07_csv/Compound.csv")

public_id <- 
compound$public_id

for(i in 50001:length(public_id)){
  cat(i, " ")
  result <- 
    tryCatch(request_foodb_compound(compound_id = public_id[i]), error = function(e) NULL)
  if(is.null(result)){
    next()
  }
  save(result, file = file.path("compound", public_id[i]), compress = "xz")
}




content <-
  read_csv("foodb_2020_04_07_csv/Content.csv")

food =
  read_csv("foodb_2020_04_07_csv/Food.csv")

food_taxonomy =
  read_csv("foodb_2020_04_07_csv/FoodTaxonomy.csv")

accession_number =
  read_csv("foodb_2020_04_07_csv/AccessionNumber.csv")





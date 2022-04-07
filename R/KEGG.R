masstools::setwd_project()
rm(list= ls())
setwd("raw_data/HMDB/")

library(metid)
library(tidyverse)

load("hmdb_database0.0.3.rda")

load("hmdb_ms1_database0.0.3.rda")

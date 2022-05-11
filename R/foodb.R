library(tidyverse)

####extract_version


####extract_metabolites
library(rvest)
request_foodb_compound_info <-
  function(url = "https://foodb.ca/compounds",
           sleep = 1) {
    result1 <-
      rvest::read_html(x = url) %>%
      rvest::html_table(fill = TRUE) %>%
      `[[`(1)
    
    result <-
      purrr::map(2:2838, function(idx) {
        cat(idx, " ")
        Sys.sleep(sleep)
        new_url <-
          paste0(url, "?page=", idx)
        x <-
          rvest::read_html(x = new_url)
        
        x <-
          x %>%
          rvest::html_table(fill = TRUE) %>%
          `[[`(1)
        x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    out <-
      rbind(result1,
            result)
    
    invisible(out)
  }


####extract_metabolites
request_foodb_compound <-
  function(url = "https://foodb.ca/compounds",
           compound_id = "FDB000004",
           return_form = c("data.frame", "list")) {
    return_form <- match.arg(return_form)
    result <-
      readLines(paste0(url, "/", compound_id, ".xml"))
    result <-
      XML::xmlTreeParse(file = result, asText = TRUE)
    result <-
      XML::xmlToList(result)
    
    names(result)
    
    if (return_form == "list") {
      result <- result
    } else{
      result <- result
    }
    invisible(result)
  }

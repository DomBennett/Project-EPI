##Libraries
library(stringr)
library(plyr)
library(RJSONIO)

taxaNames <- function(name, return.method, full.details = TRUE, limit = 10){
  # Search taxomnomic names in paleoDB
  #
  # Arguments
  #  name: string of name to be searched
  #  return.method: string specifying the return method (single, children, parents or auto)
  #  full.details: return all details or not (Boolean)
  #  limit: number of returned records (auto method only)
  #
  # Return
  #  Parsed JSON for the R envrionment
  # TODO: create dataframe of JSON object, add full detail bool, allow IDs to be used,
  #  find safe url package in R (e.g. what if there is non-url friendly characters in name
  #  given?)
  url <- "http://paleobiodb.org/data1.1/taxa/"
  if (return.method == "single") {
    url <- paste0(url, "single.json?name=", name)
  } else if (return.method == "children") {
    url <- paste0(url, "list.json?name=", name, "&rel=all_children")
  } else if (return.method == "parents") {
    url <- paste0(url, "list.json?name=", name, "&rel=all_parents")
  } else if (return.method == "auto") {
    url <- paste0(url, "auto.json?name=", name, "&limit=", limit)
  } else {
    stop("Invalid return.method specified")
  }
  res <- fromJSON(url)
  return (res)
}

# Test
taxaNames(name = "homo sapiens", return.method = 'auto')
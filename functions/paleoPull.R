##Libraries
library(stringr)
library(plyr)
library(RJSONIO)

names <- "Orycteropus"
url <- "http://paleobiodb.org/data1.1/taxa/list.json?"
names2 <- paste("name=", paste(str_replace_all(names, " ", "+"), collapse = "|", sep = ""), sep = "")
query <- paste(compact(list(url, names2)), collapse = "")
#search via API
fromJSON(query)

test.query <- "http://paleobiodb.org/data1.1/taxa/single.json?name=Dascillidae&show=attr"
fromJSON(test.query)
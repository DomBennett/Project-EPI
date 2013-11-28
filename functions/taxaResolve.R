##Libraries
library(stringr)
library(plyr)
library(RJSONIO)

taxaResolve <- function (names){
  #create query from names
  url <- "http://resolver.globalnames.org/name_resolvers.json?"
  data_source_ids <- "&data_source_ids=4" # 4 is for NCBI, only working with NCBI at the moment but you can change this.
  names2 <- paste("names=", paste(str_replace_all(names, " ", "+"), collapse = "|", sep = ""), sep = "")
  query <- paste(compact(list(url, names2, data_source_ids)), collapse = "")
  #search via API
  data <- fromJSON(query)$data
  #transform results into output
  search.name <- name.string <- canonical.form <- lineage <- 
    lineage.ids <- taxid <- match.type <- prescore <- score <- rep(NA, length(names))
  for(i in 1:length(names)){
    #print(i)
    if(length(data[[i]]) == 1){
      search.name[i] <- data[[i]][[1]]
      name.string[i] <- canonical.form[i] <- lineage[i] <- lineage.ids[i] <- taxid[i] <- 
        match.type[i] <- prescore[i] <- score[i] <- NA
    } else {
      search.name[i] <- data[[i]][[1]]
      name.string[i] <- data[[i]]$results[[1]]$name_string
      canonical.form[i] <- data[[i]]$results[[1]]$canonical_form
      lineage[i] <- data[[i]]$results[[1]]$classification_path
      lineage.ids[i] <- data[[i]]$results[[1]]$classification_path_ids
      taxid[i] <- data[[i]]$results[[1]]$taxon_id
      match.type[i] <- data[[i]]$results[[1]]$match_type
      prescore[i] <- data[[i]]$results[[1]]$prescore
      score[i] <- data[[i]]$results[[1]]$score
    }
  }
  res <- data.frame(search.name = search.name, name.string = name.string, canonical.form = canonical.form,
                    lineage = lineage, lineage.ids = lineage.ids, taxid = taxid, match.type = match.type,
                    prescore = prescore, score = score)
  return(res)
}
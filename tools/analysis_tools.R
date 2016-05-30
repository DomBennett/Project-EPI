# LIBS
library(stringdist)
library(wordcloud)
library(tm)

cateAsNum <- function(cate) {
  # Convert IUCN redlist category to number
  if(cate == 'LC') {
    return(0)
  }
  if(cate == 'NT') {
    return(1)
  }
  if(cate == 'VU') {
    return(2)
  }
  if(cate == 'EN') {
    return(3)
  }
  if(cate == 'CR') {
    return(4)
  }
  if(cate == 'EW') {
    return(5)
  }
  if(cate == 'EX') {
    return(6)
  }
  NA
}

getHbttEclgy <- function(json_lst) {
  # extract habitat and ecology descriptions from returned JSON
  hbtts <- rep(NA, length(json_lst))
  for(i in 1:length(json_lst)) {
    if(length(json_lst[[i]][['result']]) > 0) {
      res <- json_lst[[i]][['result']][[1]][['habitat']]
      if(!is.null(res)) {
        res <- gsub("<.*?>", "", res)  # remove HTML
        hbtts[i] <- res
      }
    }
  }
  names(hbtts) <- names(json_lst)
  hbtts <- hbtts[!is.na(hbtts)]
  hbtts
}

calcStrDst <- function(txts, mthd) {
  # calculate the distance between vector of texts
  cmbs <- combn(length(txts), 2)
  dsts <- rep(NA, ncol(cmbs))
  for(i in 1:ncol(cmbs)) {
    a <- txts[cmbs[ ,i]][1]
    b <- txts[cmbs[ ,i]][2]
    dsts[i] <- stringdist(a, b, method=mthd)
  }
  res <- data.frame('min'=min(dsts), 'max'=max(dsts),
                    'mean'=mean(dsts), 'sd'=sd(dsts),
                    'median'=median(dsts))
  res
}

getWrdFrq <- function(txts, min_wrd_sz=5, min_freq=1) {
  # Return the frequency of unique words in each string of vector txts
  cleanWrds <- function(txt) {
    wrds <- strsplit(txt, " ")[[1]]
    wrds <- tolower(wrds)
    wrds <- removePunctuation(wrds)
    wrds <- removeNumbers(wrds)
    wrds <- wrds[nchar(wrds) >= min_wrd_sz]
    wrds <- sub('es$', "", wrds)  # remove plurals
    wrds <- sub('s$', "", wrds)
    wrds <- unique(wrds)
    wrds
  }
  wrds <- NULL
  for(i in 1:length(txts)) {
    wrds <- c(wrds, cleanWrds(txts[i]))
  }
  wrds <- table(wrds)
  wrds <- wrds[wrds > min_freq]
  wrds
}

genTrmMtrx <- function(lfs_hes, nlfs_hes, min_wrd_sz=5, min_freq=1) {
  cleanWrds <- function(txt) {
    wrds <- strsplit(txt, " ")[[1]]
    wrds <- tolower(wrds)
    wrds <- removePunctuation(wrds)
    wrds <- removeNumbers(wrds)
    wrds <- wrds[nchar(wrds) >= min_wrd_sz]
    wrds <- sub('es$', "", wrds)  # remove plurals
    wrds <- sub('s$', "", wrds)
    wrds <- unique(wrds)
    wrds
  }
  lf_wrds <- NULL
  for(i in 1:length(lfs_hes)) {
    lf_wrds <- c(lf_wrds, cleanWrds(lfs_hes[i]))
  }
  lf_wrds <- table(lf_wrds)
  lf_wrds <- lf_wrds[lf_wrds > min_freq]
  null_wrds <- NULL
  for(i in 1:length(nlfs_hes)) {
    null_wrds <- c(null_wrds, cleanWrds(nlfs_hes[i]))
  }
  null_wrds <- table(null_wrds)
  null_wrds <- null_wrds[null_wrds > min_freq]
  # convert to matrix
  rnms <- sort(unique(c(names(null_wrds), names(lf_wrds))))
  tm <- matrix(0, ncol=2, nrow=length(rnms))
  colnames(tm) <- c('living fossil', 'null')
  rownames(tm) <- rnms
  tm[match(names(lf_wrds), rnms), 1] <- lf_wrds
  tm[match(names(null_wrds), rnms), 2] <- null_wrds
  tm
}
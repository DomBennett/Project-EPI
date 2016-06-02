# LIBS
library(stringdist)

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

justText <- function(txt) {
  # converts text to its readable form by removing punctuation and numbers
  gsub("[^a-zA-Z ]", "", txt)
}

calcStrDst <- function(txts) {
  # calculate the distance between vector of texts
  # or for a random subset of texts if a list of list of texts
  if(class(txts) == "list") {
    txts <- unlist(lapply(txts, function(x) sample(x, 1)))
  }
  cmbs <- combn(length(txts), 2)
  jw_dsts <- lv_dsts <- cosine_dsts <- rep(NA, ncol(cmbs))
  for(i in 1:ncol(cmbs)) {
    a <- justText(txts[cmbs[ ,i]][[1]])
    b <- justText(txts[cmbs[ ,i]][[2]])
    lv_dsts[i] <- stringdist(a, b, method="lv")
    cosine_dsts[i] <- stringdist(a, b, method="cosine")
    jw_dsts[i] <- stringdist(a, b, method="jw")
  }
  res <- data.frame('lv_min'=min(lv_dsts),
                    'lv_max'=max(lv_dsts),
                    'lv_mean'=mean(lv_dsts),
                    'lv_sd'=sd(lv_dsts),
                    'lv_median'=median(lv_dsts),
                    'cosine_min'=min(cosine_dsts),
                    'cosine_max'=max(cosine_dsts),
                    'cosine_mean'=mean(cosine_dsts),
                    'cosine_sd'=sd(cosine_dsts),
                    'cosine_median'=median(cosine_dsts),
                    'jw_min'=min(jw_dsts),
                    'jw_max'=max(jw_dsts),
                    'jw_mean'=mean(jw_dsts),
                    'jw_sd'=sd(jw_dsts),
                    'jw_median'=median(jw_dsts))
  res
}

getWrdFrq <- function(txts, wts=rep(1, length(txts)), min_wrd_sz=5, min_freq=5) {
  # Return the frequency of unique words in each string of vector txts
  # Use wts to determine the wt of each txt (e.g. a single clade is represented by mutliple texts)
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
  wrds <- vector("list", length=length(txts))
  for(i in 1:length(txts)) {
    wrds[[i]] <- cleanWrds(txts[i])
  }
  nreps <- unlist(lapply(wrds, function(x) length(x)))
  wts <- rep(wts, times=nreps)
  wrds <- unlist(wrds)
  tbl_wrds <- tapply(wts, wrds, sum)
  tbl_wrds <- tbl_wrds[tbl_wrds > min_freq]
  sort(tbl_wrds, decreasing=TRUE)
}
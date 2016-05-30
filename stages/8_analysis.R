# TEST IF LIVING FOSSILS SHARE MORE

# START
cat(paste0('\nStage `analysis` started at [', Sys.time(), ']\n'))

# PARAMETERS
source('parameters.R')

# FUNCTIONS
source(file.path('tools', 'analysis_tools.R'))

# DIRS
output_dir <- '8_analysis'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_dir <- "7_iucn"
output_file <- file.path("8_analysis", "res.RData")
iucn_files <- list.files(input_dir)

# LOOP THROUGH ANALYSIS GROUPS
for(iucn_file in iucn_files) {
  # INPUT
  grp <- sub("\\.RData", "", iucn_file)
  cat('    Working on [', grp, '] ....\n', sep="")
  load(file.path(input_dir, iucn_file))
  res <- matrix(ncol=6, nrow=4)
  colnames(res) <- c("Obs_mean", "Obs_sd", "Null_mean",
                     "Null_sd", "Z_score", "P_val")
  rownames(res) <- c("Categories", "Nhabitats", "Ncountries",
                     "Desciptions")
  pdf(file.path(output_dir, paste0(grp, ".pdf")))
  
  # CATEGORIES
  cat('    Testing categories ....\n')
  cates <- vector(length=length(lf_data))
  for(i in 1:length(lf_data)) {
    tmp <- vector(length=length(lf_data[[i]]))
    for(j in 1:length(lf_data[[i]])) {
        tmp[j] <- cateAsNum(lf_data[[i]][[j]][['cate']])
    }
    cates[i] <- mean(tmp, na.rm=TRUE)
  }
  null <- vector(length=length(null_cate))
  for(i in 1:length(null_cate)) {
    tmp <- unlist(lapply(null_cate[[i]],
                         function(x) cateAsNum(x)))
    null[i] <- mean(tmp, na.rm=TRUE)
  }
  obs_mean <- mean(cates, na.rm=TRUE)
  obs_sd <- sd(cates, na.rm=TRUE)
  null_mean <- mean(null, na.rm=TRUE)
  null_sd <- sd(null, na.rm=TRUE)
  p_val <- sum(null <= obs_mean)/length(null)
  z_score <- (obs_mean - null_mean)/null_sd
    sum(null <= obs_mean)/length(null)
  hist(null, main=paste0('P = ', signif(p_val, 3)),
       xlab="Categories")
  abline(v=obs_mean, col='red')
  res[1, ] <- c(obs_mean, obs_sd, null_mean,
                null_sd, z_score, p_val)
  cat("    Done.\n")
  
  # HABITATS
  cat('    Testing habitats ....\n')
  nhbbts <- vector(length=length(lf_data))
  for(i in 1:length(lf_data)) {
    tmp <- vector(length=length(lf_data[[i]]))
    for(j in 1:length(lf_data[[i]])) {
      tmp[j] <- lf_data[[i]][[j]][['nhbbts']]
    }
    nhbbts[i] <- mean(tmp, na.rm=TRUE)
  }
  null <- vector(length=length(null_nhbbts))
  for(i in 1:length(null_nhbbts)) {
    tmp <- unlist(lapply(null_nhbbts[[i]],
                         function(x) x))
    null[i] <- mean(tmp, na.rm=TRUE)
  }
  obs_mean <- mean(nhbbts, na.rm=TRUE)
  obs_sd <- sd(nhbbts, na.rm=TRUE)
  null_mean <- mean(null, na.rm=TRUE)
  null_sd <- sd(null, na.rm=TRUE)
  p_val <- sum(null <= obs_mean)/length(null)
  z_score <- (obs_mean - null_mean)/null_sd
  sum(null <= obs_mean)/length(null)
  hist(null, main=paste0('P = ', signif(p_val, 3)),
       xlab="Habitats")
  abline(v=obs_mean, col='red')
  res[2, ] <- c(obs_mean, obs_sd, null_mean,
                null_sd, z_score, p_val)
  cat("    Done.\n")
  
  # HABITATS
  cat('    Testing countries ....\n')
  ncntrs <- vector(length=length(lf_data))
  for(i in 1:length(lf_data)) {
    tmp <- vector(length=length(lf_data[[i]]))
    for(j in 1:length(lf_data[[i]])) {
      tmp[j] <- lf_data[[i]][[j]][['ncntrs']]
    }
    ncntrs[i] <- mean(tmp, na.rm=TRUE)
  }
  null <- vector(length=length(null_ncntrs))
  for(i in 1:length(null_ncntrs)) {
    tmp <- unlist(lapply(null_ncntrs[[i]],
                         function(x) x))
    null[i] <- mean(tmp, na.rm=TRUE)
  }
  obs_mean <- mean(ncntrs, na.rm=TRUE)
  obs_sd <- sd(ncntrs, na.rm=TRUE)
  null_mean <- mean(null, na.rm=TRUE)
  null_sd <- sd(null, na.rm=TRUE)
  p_val <- sum(null <= obs_mean)/length(null)
  z_score <- (obs_mean - null_mean)/null_sd
  sum(null <= obs_mean)/length(null)
  hist(null, main=paste0('P = ', signif(p_val, 3)),
       xlab="Countries")
  abline(v=obs_mean, col='red')
  res[3, ] <- c(obs_mean, obs_sd, null_mean,
                null_sd, z_score, p_val)
  cat("    Done.\n")
  
  # DESCRIPTION DIFFERENCE
  cat('    Testing description difference ....\n')
  lf_dst <- matrix(nrow=nitrtns, ncol=5)
  ns <- vector(length=length(nitrtns))
  for(itrtn in 1:nitrtns) {
    txts <- vector(length=length(lf_data))
    for(i in 1:length(lf_data)) {
      j <- sample(1:length(lf_data[[i]]), 1)
      txt <- lf_data[[i]][[j]][['dscrptn']]
      if(!is.null(txt)) {
        txts[j] <- gsub("<.*?>", "", txt)  # remove HTML
      }
    }
    bool <- txts != "FALSE"
    ns[itrtn] <- sum(!bool)
    tmp <- calcStrDst(txts[bool], mthd='cosine')
    lf_dst[itrtn, ] <- as.numeric(tmp)
  }
  colnames(lf_dst) <- colnames(tmp)
  null_dsts <- matrix(nrow=nitrtns, ncol=5)
  for(i in 1:length(null_dscrptn)) {
    txts <- as.character(unlist(null_dscrptn[[i]]))
    txts <- sample(txts, ns[i], replace=TRUE)
    txts <- gsub("<.*?>", "", txts)
    tmp <- calcStrDst(txts, mthd="cosine")
    null_dsts[i, ] <- as.numeric(tmp)
  }
  colnames(null_dsts) <- colnames(tmp)
  obs_mean <- mean(lf_dst[,'median'], na.rm=TRUE)
  obs_sd <- sd(lf_dst[,'median'], na.rm=TRUE)
  null_mean <- mean(null_dsts[, 'median'], na.rm=TRUE)
  null_sd <- sd(null_dsts[, 'median'], na.rm=TRUE)
  p_val <- sum(null_dsts[, 'median'] <= obs_mean)/nrow(null_dsts)
  z_score <- (obs_mean - null_mean)/null_sd
  hist(null_dsts, main=paste0('P = ', signif(p_val, 3)),
       xlab="Description")
  abline(v=obs_mean, col='red')
  res[4, ] <- c(obs_mean, obs_sd, null_mean,
                null_sd, z_score, p_val)
  cat("    Done.\n")
  dev.off()
  
  # OUTPUT
  save(res, file=file.path(output_dir, iucn_file))
}

# END
cat(paste0('\nStage `analysis` finished at [', Sys.time(), ']\n'))
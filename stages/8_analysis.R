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
    txts <- sample(txts, ns[i], replace=TRUE)  # TODO: should this really be replace?
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
  
  # WORD FREQUENCIES
  wts <- txts <- NULL
  for(i in 1:length(lf_data)) {
    cc <- 0
    for(j in 1:length(lf_data[[i]])) {
      txt <- lf_data[[i]][[j]][['dscrptn']]
      if(!is.null(txt)) {
        txts <- c(txts, gsub("<.*?>", "", txt))
        cc <- cc + 1
      }
    }
    wts <- c(wts, rep(1/(cc+1), cc))
  }
  obs_frq <- getWrdFrq(txts, wts, min_freq=2)
  null_frqs <- vector("list", length=nitrtns)
  for(i in 1:nitrtns) {
    txts <- null_dscrptn[[i]]
    txts <- gsub("<.*?>", "", txts)
    null_frqs[[i]] <- getWrdFrq(txts, min_freq=2)
  }
  # how many more times does a term appear in observed?
  frq_res <- data.frame(wrd=NA, obs_freq=NA, exp_freq=NA,
                        p_val=NA, z_score=NA)
  for(i in 1:length(obs_frq)) {
    wrd <- names(obs_frq)[i]
    nd <- unlist(lapply(null_frqs, function(x) {
      if(wrd %in% names(x)) {
        res <- x[[wrd]]
      } else {
        res <- 0
      }
      res
    }))
    frq_res[i, 'wrd'] <- wrd
    frq_res[i, 'obs_freq'] <- obs_frq[[wrd]]
    frq_res[i, 'exp_freq'] <- mean(nd)
    frq_res[i, 'z_score'] <- (obs_frq[[wrd]] - mean(nd))/sd(nd)
    frq_res[i, 'p_val'] <- sum(nd >= obs_frq[[wrd]])/length(nd)
  }
  plot_res <- frq_res[frq_res$p_val < 0.05, ]
  plot_res <- plot_res[plot_res$z_score < 5, ]
  pdf(file.path(output_dir, paste0(grp, "_wordcloud.pdf")), w=14, h=14)
  wordcloud(plot_res$wrd, plot_res$z_score,
            colors=brewer.pal(8, 'Dark2'), max.words=45)
  dev.off()
  
  # OUTPUT
  save(res, frq_res, file=file.path(output_dir, iucn_file))
}

# END
cat(paste0('\nStage `analysis` finished at [', Sys.time(), ']\n'))
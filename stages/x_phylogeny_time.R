library(treeman)

# DIRS
data_dir <- file.path("0_data", "raw")
nt_dir <- "x_ncbi_taxonomy"

load(file.path(nt_dir, "res.Rdata"))

# IDENTIFY MATCHING NODES

# MATCH TIPS TO NMS
tps <- gsub("_", " ", tree['tips'])
for(i in 1:length(nms_obj)) {
  bool <- nms_obj[[i]] %in% tps
  if(any(bool)) {
    node_obj[[i]][['tip_i']] <- which(tps == nms_obj[[i]][bool][1])
  }
}

# MATCH NODES OF INTEREST TO NODES IN TREE
txids <- names(nms_obj)
tree_kids <- getNdsKids(tree, ids=tree['nds'])
for(i in which(!ignore_bool)) {
  nd_kids <- node_obj[[i]][['kids']]
  if(is.null(nd_kids)) {
    next
  }
  if(nd_kids[1] == "none") {
    tip_i <- node_obj[[i]][['tip_i']]
    if(is.null(tip_i)) {
      next
    }
    node_obj[[i]][["nid"]] <- tree['tips'][tip_i]
  } else {
    nd_tips <- vector(length=length(nd_kids))
    for(j in 1:length(nd_kids)) {
      tip_i <- node_obj[[which(txids == nd_kids[j])]][['tip_i']]
      if(!is.null(tip_i)) {
        nd_kids[j] <- tree['tips'][tip_i]
      }
    }
    mtch_scrs <- vector(length=length(tree_kids))
    for(j in 1:length(tree_kids)) {
      mtch_scrs[j] <- (sum(nd_kids %in% tree_kids[[j]])/length(nd_kids)) +
        (sum(tree_kids[[j]] %in% nd_kids)/length(tree_kids[[j]]))
    }
    bst_mtch <- which.max(mtch_scrs)[1]
    if(length(bst_mtch) > 0) {
      node_obj[[i]][["nid"]] <- tree['nds'][bst_mtch]
    }
  }
}

# FIND PD, ED, SPAN AND AGE FOR EVERY NODE
ed_vals <- calcFrPrp(tree, tids=tree['tips'])
for(i in which(!ignore_bool)) {
  nid <- node_obj[[i]][['nid']]
  if(is.null(nid)) {
    next
  }
  age <- getNdAge(tree, id=nid)
  node_obj[[i]][['age']] <- age
  if(age < 50) {
    ignore_bool[i] <- TRUE
    next
  }
  kids <- getNdKids(tree, nid)
  ed <- mean(ed_vals[which(tree['tips'] %in% kids)])
  node_obj[[i]][['ed']] <- ed
  spn <- getNdSlt(tree, slt_nm="spn", id=nid)
  node_obj[[i]][['spn']] <- spn
  pd <- getNdSlt(tree, slt_nm="pd", id=nid)
  node_obj[[i]][['pd']] <- pd
}

# TOP-100
spns <- unlist(lapply(node_obj, function(x) x[['spn']]))
cns <- unlist(lapply(node_obj, function(x) suppressWarnings(min(x[['cntrst_n']]))))
nms <- unlist(lapply(nms_obj, function(x) x[['scientific name']]))
bool <- unlist(lapply(node_obj, function(x) !is.null(x[['spn']])))
epis <- log(cns[bool]) * log(spns)
ordrd <- order(epis)[1:100]
for(i in ordrd) {
  nm <- nms[i]
  epi <- epis[i]
  spcr <- paste0(rep(' ', 33 - nchar(nm)), collapse="")
  cat(nm, spcr, signif(epi, 3), "\n")
}

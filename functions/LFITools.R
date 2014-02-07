## No copyright, no warranty
## Dominic John Bennett
## Functions for calculating LFI
## 07/02/2014

## Dependencies
require (ape)
require (geiger)

## Functions
listClades <- function (phylo) {
  # List all descendants of each node in phylogeny.
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #
  # Returns:
  #  list of list of clades and node number
  nodes <- 1:(phylo$Nnode + length (phylo$tip.label))
  nodes <- nodes[-(length (phylo$tip.label) + 1)]
  clades <- list()
  for (node in nodes) {
    clades <- c (clades, list (tips(phylo, node)))
  }
  return (list (clades, nodes))
}

matchClades <- function (qclades, sclades, exact = FALSE) {
  # Phylogenies may have different numbers of species, as such node numbers will not be
  #  directly comparable. This function matchs nodes between phylogenies using clade names.
  #
  # Args:
  #  qclades: query clades, a list of clades, the clades to be mapped to the subject clades
  #  sclades: subject clades
  #  exact: if true, nodes will be mapped to one another if they have all of the same
  #   descendants. If false, nodes will be mapped if descedants are missing to allow for 
  #   missing taxa (i.e. smallest node with all the same descendants).
  #
  # Returns:
  #  vector of query node numbers
  if (exact) {
    return (match (qclades, sclades))
  } else {
    subFunction <- function (qclade, clades) {
      scores <- unlist (lapply (clades, function (x) sum (qclade %in% x)))
      if (sum (scores) == 0) {
        return (NA)
      } else {
        res <- clades[scores == max (scores)]
        len.res <- sapply (res, length)
        res <- res[len.res == min (len.res)]
        return (match (res, clades))
      }
    }
    return (unlist (lapply (qclades, subFunction, clades = sclades)))
  }
}

addOutgroup <- function (phylo, outgroup.factor = 100) {
  # Add an outgroup to phylogeny
  #  (https://stat.ethz.ch/pipermail/r-sig-phylo/2012-November/002417.html)
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  outgroup.factor: determines length of branch connecting outgroup to the rest of phylo,
  #   large values produce smaller branches (default 100)
  #
  # Return:
  #  phylogeny (ape class)
  rtt.dist <- mean (diag (vcv.phylo (phylo)))
  dist <- rtt.dist/outgroup.factor
  edge.matrix <- matrix (c (3,2,3,1), 2, 2, byrow = TRUE)
  tip <- list (edge = edge.matrix, tip.label = c ("outgroup", "clade.to.be"),
              edge.length = c (dist, rtt.dist + dist), Nnode = 1)
  class (tip) <- "phylo"
  return (bind.tree (tip, phylo, where = 2))
}

parsimonyReconstruction <- function (chars, phylo, missing.char = "?",
                                        order.numeric = TRUE) {
  # Return upper and lower character states for all nodes in phylogeny for a list of
  #  characters using parsimony. Where missing data occur, the phylogeny is reduced by
  #  dropping tips. The ancestral character is determined as the most primitive state
  #  for ordered characters and randomly for unordered characters.
  #
  # Args:
  #  chars: matrix with species names as rows and character states as columns
  #  phylo: phylogeny (ape class)
  #  missing.char: symbol for missing data (default "?")
  #  order.numeric: if true, assumes all numeric characters are ordered
  #
  # Returns:
  #  list of node states and reduced phylogeny for each character (reconstruction.obj)
  #if (!is.ultrametric (phylo)) {
  #  stop("Phylogeny must be ultrametric (i.e. w/o polytomies)")
  #}
  if (all (rownames (chars) %in% phylo$tip.label)) {
    stop("Chars must be a matrix w/ its rownames matching tips in phylogeny")
  }
  reconstruction.obj <- list ()
  for (i in 1:ncol (chars)) {
    temp.chars <- chars[phylo$tip.label,i]
    tips.to.drop <- names (temp.chars)[temp.chars == missing.char]
    temp.chars <- temp.chars[temp.chars != missing.char]
    if (length (tips.to.drop) > 0) {
      reduced.tree <- drop.tip (phylo, tips.to.drop)
    } else {
      reduced.tree <- phylo
    }
    # Convert chars to numeric and use primitive state for outgroup
    non.numeric.chars <- FALSE
    if (order.numeric) {
      temp.chars.names <- names (temp.chars)
      temp.chars <- tryCatch (as.numeric (temp.chars),
                             warning = function (c) temp.chars)
      if(is.numeric (temp.chars)) {
        names (temp.chars) <- temp.chars.names
        reduced.tree <- addOutgroup (reduced.tree)
        temp.chars <- c (temp.chars, min (temp.chars))
        # min is the most primitive state
        names (temp.chars)[length (temp.chars)] <- "outgroup"
      } else {
        non.numeric.chars <- TRUE
      }
    }
    # Else, use random character for outgroup
    if (!order.numeric | non.numeric.chars) {
      reduced.tree <- addOutgroup (reduced.tree)
      temp.chars <- c (temp.chars, sample (temp.chars, 1))
      names (temp.chars)[length (temp.chars)] <- "outgroup"
    }
    if (is.rooted (reduced.tree)) {
      reduced.tree <- unroot (reduced.tree)
    }
    # calculate for nodes using parsimony
    res <- MPR (temp.chars[reduced.tree$tip.label], reduced.tree, "outgroup")
    # add tip node states
    res <- rbind (matrix (rep (temp.chars[reduced.tree$tip.label], each = 2), ncol = 2,
                        byrow = TRUE), res)
    res <- list (res, reduced.tree)
    reconstruction.obj <- c (reconstruction.obj, list (res))
  }
  return (reconstruction.obj)
}

calcBranchChanges <- function (phylo, reconstruction.obj, as.mean = TRUE, ...) {
  # Count the number of changes that have occured along each branch in phylogeny given
  #  a reconstruction object.
  #
  # Args:
  #  phylo: full phylogeny (ape class)
  #  reconstruction.obj: a list for multiple characters containing a matrix of upper and lower
  #   estimates of node states and a reduced phylogeny (i.e. return from
  #   parsominyReconstruction)
  #  as.mean: return a vector for the mean change per branch, or a list for each change
  #   per character per branch (default true)
  #
  # Returns:
  #  vector or list
  clades <- listClades (phylo)[[1]]
  nodes <- listClades (phylo)[[2]]
  res <- list()
  for (i in 1:length (reconstruction.obj)) {
    node.states <- reconstruction.obj[[i]][[1]]
    reduced.tree <- reconstruction.obj[[i]][[2]]
    part.res <- matrix (rep (0, length (clades) * 2), nrow = 2)
    rownames (part.res) <- c("Change", "Node.presence")
    colnames (part.res) <- nodes
    temp.clades <- listClades (reduced.tree)[[1]]
    temp.nodes <- listClades (reduced.tree)[[2]]
    # find shared nodes between phylo and reduced tree
    clade.index <- matchClades (temp.clades, clades, ...)
    edge.store <- edge.score <- vector ()
    for (j in 1:length (clade.index)) {
      if (is.na (clade.index[j])) {
        next
      } else {
        edge.store <- c (edge.store, which (reduced.tree$edge[ ,2] == temp.nodes[j]))
        edge <- reduced.tree$edge[edge.store[length (edge.store)], ]
        start <- character[edge[1], ]
        end <- character[edge[2], ]
        if (is.numeric(start)) {
          # for ordered states
          change <- abs (sum (start) - sum (end))/2 
        } else {
          # for unordered states
          change <- sum (start != end)/2
        }
        edge.score <- c (edge.score, change)
        part.res[1, clade.index[j]] <- change
        part.res[2, clade.index[j]] <- 1
      }
    }
    res <- c (res, list (part.res))
  }
  if (as.mean) {
    tot.changes <- rowSums(matrix (unlist (lapply (res, function (x) x [1, ])),
                           nrow = length (clades), byrow = TRUE))
    tot.nchars <- rowSums(matrix (unlist (lapply (res, function (x) x [2, ])),
                                   nrow = length (clades), byrow = TRUE))
    res <- tot.changes/nchars
    names (res) <- nodes
  }
  return (res)
}

plotBranchChanges <- function (phylo, node.states, changes) {
  # Plot phylogeny with character labels on tips, edges and nodes
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  node.states: matrix, upper and lower bound estimate for each node
  #  changes: changes for each edge
  #
  # Returns:
  #  None
  plot(phylo, show.tip.label = FALSE, no.margin = TRUE)
  nodelabels(paste("[", node.states[-(1:length(phylo$tip.label)), 1], ",",
                   node.states[-(1:length(phylo$tip.label)), 2], "]", sep = ""))
  tiplabels(node.states[1:length(phylo$tip.label),1], adj = -2)
  edgelabels(text = changes, edge = names(changes), adj = c(0, 1))
}
  

calcLFI <- function (phylo, branch.changes) {
  # Calculate a living fossil index based on branch changes
  return (NA)
}

plotLFI <- function (phylo, lfi.res) {
  # Plot the results of LFI
  return (NA)
}

## 17/06/2014
## D.J. Bennett
## Do clades rise and fall?
# Here I use an equal rates markov model to add species to
# a tree. I then record the success of each clade at set
# time intervals. Since I have not added any extinction
# I want to see if each clade shows a signal of rising.
# Future: add extinction, time steps are linked to the
# total number of species - have a random number of new tips,
# record tree growth in a video

library (ape)
library (MoreTreeTools)
library (plyr)

countChildren <- function (tree) {
  # Count the number of children for every node
  .count <- function (node.label) {
    node <- which (tree$node.label == node.label)
    node <- length (tree$tip.label) + node
    length (getChildren (tree, node))
  }
  # all internal nodes
  node.labels <- tree$node.label
  res <- mdply (.data = data.frame (node.label = node.labels),
                .fun = .count)
  colnames (res) <- c ('node', 'n.children')
  res
}

reformat <- function (clade.performance) {
  # take list of list and convert to a dataframe
  .getTime <- function (i, node) {
    # get success for node at a time point
    data <- clade.performance[[i]]
    if (any (data$node == node)) {
      return (data[data$node == node, 2])
    }
    0
  }
  .getNode <- function (node) {
    mdply (.data = data.frame (
      i = 1:length (clade.performance)),
      .fun = .getTime, node)[ ,2]
  }
  .addNode <- function (node) {
    # add success for node at all time points
    # for a res dataframe
    node <- as.character (node)
    node.success <- .getNode (node)
    res[node] <- node.success
    res <<- res
  }
  # get nodes across times
  nodes <- unique (unlist (llply (.data = clade.performance,
                  .fun = function (x) as.vector(x$node))))
  # build res dataframe by adding first results
  res <- data.frame (.getNode (nodes[1]))
  colnames (res) <- nodes[1]
  nodes <- data.frame (node = nodes[-1])
  # add to res
  m_ply (.data = nodes, .fun = .addNode)
  res
}

plotSuccess <- function (res, section = 'all',
                         thinning = 'all') {
  # Plot as separate lines for the success of each
  # clade through time
  # section = 'all' -- for all clads
  # else do for bottom number
  # thinning allows only a proportion of the
  # clade results to be plotted
  if (section != 'all') {
    top <- ncol (res)
    bottom <- top - section
    res <- res[bottom:top]
  }
  if (thinning != 'all') {
    # choose those to plot by biasing towards the
    # older ones -- makes it easier to see
    # too many young ones
    random.is <- sample (1: ncol (res),
                         thinning, prob = 1/1: ncol (res))
    res <- res[ ,random.is]
  }
  n.clades <- ncol (res)
  # re-use the first 3 rainbow colours
  # red, green and blue
  nreps <- ceiling (n.clades/3)
  cols <- rep (rainbow (3, alpha = 0.7), nreps)
  plot (res[,1], col = NULL, ylab = 'N', xlab = 'Time',
        main = 'Clade success through time')
  for (i in 1:n.clades) {
    non.zero <- res [,i] != 0
    y <- res[non.zero, i]
    x <- which (non.zero)
    lines (y = y, x = x, col = cols[i], pch = 19)
  }
}

growTree <- function (iterations, birth = 0.6, death = 0.4) {
  add <- function () {
    # an equal rates markov model:
    # choose a species at random
    # simulate it speciating by adding 1 to all tip edges
    # adding a new tip at a new node in the random species' connecting
    # edge at 1 time step ago
    # choose a species at random
    # grow tip.edges by 1
    tip.edges <- which (tree$edge[, 2] %in%
                          1:length (tree$tip.label))
    tree$edge.length[tip.edges] <-
      tree$edge.length[tip.edges] + 1
    #probs <- tree$edge.length[tip.edges] # weight by branch length?
    random.edge <- sample (x = tip.edges, size = 1)
    new.node.label <- paste0 ('n', max.node + 1)
    new.tip.label <- paste0 ('t', length (tree$tip.label) + 1)
    # new node.age is always 1 -- 1 time step ago.
    node.age <- 1
    # add new tip at random edge
    tree <<- addTip (tree = tree, edge = random.edge,
                     tip.name = new.tip.label,
                     node.age = node.age,
                     node.label = new.node.label)
    max.node <<- max.node + 1
  }
  drop <- function () {
    to.drop <- sample (tree$tip.label, 1)
    tree <<- drop.tip (tree, to.drop)
  }
  run <- function (iteration) {
    add.bool <- sample (c (TRUE, FALSE), size = 1,
                        prob = c (birth, death))
    if (add.bool) {
      add ()
    } else if (length (tree$tip.label) > 2) {
      drop ()
    } else {
      print ('Drop species -- but too few species in tree. Skipping.')
    }
  }
  m_ply (.data = (iteration = 1:iterations), .fun = run)
  tree
}

# parameters
time.steps <- 10000
intervals <- 10
max.node <- 1 # global counter of nodes
# seed tree of two species
tree <- stree (2)
# add lengths and labels
tree$edge.length <- c (1,1)
# N.B. ape assumes node labels only correspond to internal nodes
#  MoreTreeTools assumes they correspond to all nodes
tree$node.label <- 'n1'
iterations <- time.steps/intervals
clade.performance <- list ()
for (i in 1:iterations) {
  print (i)
  tree <- growTree (intervals, birth = 3, death = 1)
  current.success <- countChildren (tree)
  clade.performance <- c (clade.performance,
                          list (current.success))
}
res <- reformat (clade.performance)
plotSuccess (res)
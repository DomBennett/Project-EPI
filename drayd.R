## 14/01/2013
## D.J. Bennett
## Importing morpho and phylo data from dryad

## Dirs
input.dir <- "0_data/dryad"

## Libs
require(ape)

## Ladybird data: http://datadryad.org/resource/doi:10.5061/dryad.dc1r2
# charset morph = 2254-2365
data <- read.nexus.data(file.path(input.dir, "ladybird_matrix.nex"))
data <- lapply(data, function(x) x[2254:2365])
tree <- read.nexus(file.path(input.dir, "ladybird_tree.nex"))[[1]]
tree$tip.label %in% names(data)
plot(tree, "unrooted")
tree
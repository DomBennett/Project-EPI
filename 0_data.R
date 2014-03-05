## No copyright, no warranty
## Dominic John Bennett
## 04/03/2014

## Libraries
require (ape)

## Functions
matchPhyChar <- function (phylo, chars) {
  chars <- chars[rownames (chars) %in% phylo$tip.label, ]
  if (!any (phylo$tip.label %in% rownames (chars))) {
    phylo <- drop.tip (phylo,
                       tip = phylo$tip.label[!phylo$tip.label %in% rownames (chars)])
  }
  if (!all (c (phylo$tip.label %in% rownames (chars),
               rownames (chars) %in% phylo$tip.label))) {
    stop ("All labels do not match between data and phylogeny.")
  }
  return (list (phylo = phylo, chars = chars))
}


## Directories
output.dir <- "0_data"
input.dir <- file.path(output.dir, "raw")

## Ladybird
data <- read.nexus.data (file.path (input.dir, "ladybird_matrix.nex"))
data <- lapply (data, function (x) x[2254:2365])
chars <- matrix (unlist (data), nrow = length (data), byrow = TRUE)
rownames (chars) <- names (data)
rownames (chars) <- sub (" ", "_", rownames (chars))
chars[chars == "?"] <- rep (NA, sum (chars == "?"))
phylo <- read.nexus (file.path (input.dir, "ladybird_tree.nex"))[[1]]
data <- matchPhyChar (phylo, chars)
save (data, file = file.path (output.dir, "ladybird.RData"))

## Mammal (http://esapubs.org/archive/ecol/e090/184/)
data <- read.delim (file.path (input.dir, "panTHERIA.txt"), na.strings = -999,
                    stringsAsFactors = FALSE)
chars <- data[ ,- c (1:5,36:55)]
for (i in 1:ncol (chars)) {
  temp.chars <- chars[!is.na (chars[ , i]),i]
  if (length (unique (temp.chars)) > 10) {
    chars[!is.na (chars[ , i]),i] <- cut(temp.chars, 10)
  }
}
rownames (chars) <- data[ ,"MSW93_Binomial"]
rownames (chars) <- sub (" ", "_", rownames (chars))
phylo <- read.tree (file.path (input.dir, "bininda.txt"))
data <- matchPhyChar (phylo, chars)
save (data, file = file.path (output.dir, "mammal.RData"))

source(file.path('tools', 'chng_tools.R'))

library(ape)
tree <- read.tree(text="((((t7:0.09090909091,t1:0.09090909091):0.1818181818,(t8:0.09090909091,t10:0.09090909091):0.1818181818):0.4545454545,(t11:0.3636363636,((t9:0.09090909091,t12:0.09090909091):0.1818181818,(t5:0.09090909091,t3:0.09090909091):0.1818181818):0.09090909091):0.3636363636):0.2727272727,((t6:0.09090909091,t4:0.09090909091):0.09090909091,t2:0.1818181818):0.8181818182);")

plot(tree);edgelabels()


chars <- c(1, 0, 0, NA, NA, 2, 2, 4, 4, 0, 0, NA)
names(chars) <- tree$tip.label

rdcd_tree <- drop.tip(tree, tip=names(chars)[is.na(chars)])

rdcd_tree <- addOutgroup(rdcd_tree)
rdcd_tree <- unroot(rdcd_tree)
chars <- chars[!is.na(chars)]
chars <- c(chars, 0)
names(chars)[length(chars)] <- 'outgroup'
prmsmy_res <- MPR(x=chars, phy=rdcd_tree, outgroup="outgroup")
plot(rdcd_tree, show.tip.label=FALSE)
tiplabels(chars[rdcd_tree$tip.label], adj = -2)
nodelabels(paste("[", prmsmy_res[, 1], ",", prmsmy_res[, 2], "]", sep = ""))


plot(tree, show.tip.label=FALSE, edge.width=4)

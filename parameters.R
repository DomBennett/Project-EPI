# SET PARAMETERS
parallel <- TRUE
if(parallel) {
  library(doMC)
  ncps <- detectCores()
  registerDoMC(cores=ncps)
}
stdy_grp <- 'plants'

# NCBI TAXONOMY
# Division Dump Info:
# 0	|	BCT	|	Bacteria	|		|
# 1	|	INV	|	Invertebrates	|		|
# 2	|	MAM	|	Mammals	|		|
# 3	|	PHG	|	Phages	|		|
# 4	|	PLN	|	Plants and Fungi	|		|
# 5	|	PRI	|	Primates	|		|
# 6	|	ROD	|	Rodents	|		|
# 7	|	SYN	|	Synthetic and Chimeric	|		|
# 8	|	UNA	|	Unassigned	|	No species nodes should inherit this division assignment	|
# 9	|	VRL	|	Viruses	|		|
# 10	|	VRT	|	Vertebrates	|		|
# 11	|	ENV	|	Environmental samples	|	Anonymous sequences cloned directly from the environment	|
division_codes <- c(1, 2, 4, 5, 6, 10)
contrst_n_min <- 100
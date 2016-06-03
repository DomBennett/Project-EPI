# SET PARAMETERS

# PARALLEL SETTINGS
parallel <- TRUE  # Won't work for Windows
if(parallel) {
  library(doMC)
  ncps <- detectCores()
  registerDoMC(cores=ncps)
}

# DIVISION SELECTION
# division.dmp (downloaded: 16/05/2016)
# 0	 |	BCT	|	Bacteria
# 1	 |	INV	|	Invertebrates
# 2	 |	MAM	|	Mammals
# 3	 |	PHG	|	Phages
# 4	 |	PLN	|	Plants and Fungi
# 5	 |	PRI	|	Primates
# 6	 |	ROD	|	Rodents
# 7	 |	SYN	|	Synthetic and Chimeric
# 8	 |	UNA	|	Unassigned
# 9	 |	VRL	|	Viruses
# 10 |	VRT	|	Vertebrates
# 11 |	ENV	|	Environmental samples
division_codes <- c(1, 2, 4, 5, 6, 10)

# ANALYSIS GROUPS
# which taxonomic groups will be grouped for analysis?
anlyss_grps <- list("mammals"="40674",  # must match tree name in 0_data/trees
                    "birds"="8782",     # must match tree name in 0_data/trees
                    "bony_fish"="186623",
                    "saurians"="32561",
                    "amphibia"="8292",
                    "vertebrates"="7742")

# ANALYSIS PARAMETERS
nitrtns <- 999  # number of iterations for permuation test
nlfs <- 100     # number of clades at which to cutoff

# IGNORE PATTERNS
# Ignore nodes that are not true "natural" biological entities
# Ignore extinct taxa
ignore_pttrns <- c("unclassified",       # unclassified biological entities
                   "unassigned",         # unassigned biological entities
                   "\\sx\\s",            # species crosses
                   "incertae sedis",     # uncertain taxonomic group
                   "toxodon",            # extinct
                   "macrauchenia",       # extinct
                   "brachylophosaurus",  # extinct
                   "tyrannosaurus")      # extinct


# LIVING FOSSIL SETTINGS
# These determine what node is chosen as a living fossil candidate (`cnddts`).
# Not doing this will run the pipeline on all nodes, 100s to 1,000,0000s.
max_cntrst_n <- 0.01  # maximum contrasted N for a living fossil
min_prnt_n <- 1000    # minimum number of descendant species in direct parent
max_nsstrs <- 100     # maximum number of sisters, if node is polytomous
min_age <- 50         # minimum age of a node in MY
# SET PARAMETERS

# PARALLEL SETTINGS
parallel <- TRUE  # Won't work for Windows
if(parallel) {
  library(doMC)
  ncps <- detectCores() - 2
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
ignr_codes <- c(0, 3, 7, 8, 9, 11)

# ANALYSIS GROUPS
# which taxonomic groups will be grouped for analysis?
anlyss_grps <- list("mammals"="40674",  # must match tree name in 0_data/trees
                    "birds"="8782",     # must match tree name in 0_data/trees
                    "bony_fish"="186623",
                    "lepidosaurs"="8504",
                    "amphibia"="8292",
                    "vertebrates"="7742")

# URL SEARCH (timetree, ncbi, wikipedia)
wt <- c(1, 3, 10, 60, 120, 600, 3600, 7200)  # waiting times for attempts to access a URL

# TIMETREE
max_tt <- 10  # maximum number of searches if multiple sisters, higher more accurate but slower
try_again <- FALSE  # saves time split dates to prevent excessive searches,
 # if max_tt is low there is a chance the incorrect tmsplt could have been saved,
 # in which case for added accuracy (but more time) increase max_tt and set this to TRUE

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
max_cntrst_n <- 0.1  # maximum contrasted N for a living fossil
min_prnt_n <- 500    # minimum number of descendant species in direct parent
max_nsstrs <- 100     # maximum number of sisters, if node is polytomous
min_age <- 50         # minimum age of a node in MY
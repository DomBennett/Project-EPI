# DOWNLOAD IUCN DATA ON LIVING FOSSILS AND NULL SPECIES

# START
cat(paste0('\nStage `iucn download` started at [', Sys.time(), ']\n'))

# FUNCTIONS
library(gridExtra)
source(file.path('tools', 'i_tools.R'))
source(file.path('tools', 'iucn_dwnld_tools.R'))
source(file.path('tools', 'node_obj_tools.R'))
source(file.path('tools', 'synonym_tools.R'))

# PARAMETERS
source('parameters.R')
token <- getToken()
min_nchar <- 3  # minimum number of characters in a meaningful word

# DIRS
output_dir <- '9_iucn_dwnld'
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}
input_file <- file.path("8_epi", "res.RData")

# INPUT
load(input_file)
epi <- epi[!duplicated(epi$txid), ]

# LOOP FOR WORDS
wrds <- vector("list", length=nrow(epi))
cat('    Searching IUCN for descriptions ....')
for(i in 1:nrow(epi)) {
  nms <- getKidNms(epi[i,'txid'])
  res <- NULL
  for(nm in nms) {
    nrrtv <- getIUCNNrrtv(nm, token)
    if(class(nrrtv) == "list" && length(nrrtv[['result']]) > 0 &&
       !is.null(nrrtv[['result']][[1]][['habitat']])) {
      nrrtv <- nrrtv[['result']][[1]][['habitat']]
      nrrtv <- gsub("<.*?>", "", nrrtv)  # remove html tags
      res <- c(res, cleanWrds(nrrtv))
    }
  }
  if(length(res) > 0) {
    wrds[[i]] <- unique(res)
  }
}
wwrds <- which(sapply(wrds, length) > 0)
cat("Done, found data for [", length(wwrds), "/", nrow(epi),
    "].\n", sep="")

# TEST REFUGIUM
pstvs <- c("^refug", "^stable$", "^fossorial",
           "^underground$", "^static", "^dig")
ngtvs <- c("^dynamic", "^changing$", "^unstable$")
wrd_counts <- senseScore(slt_nm="refuge_scr", pstvs, ngtvs)
if(!all(unique(epi$refuge_scr) %in% c(0, -1, 1, NA))) {
  stop('Unexpected result in sense score')
}
# loop up example texts
# getExampleTextFrmWrd(epi, pttrn="^cave$")
# getExampleTextFrmI(epi, sample(which(epi$refuge_scr > 0), 1))
# getExampleTextFrmI(epi, sample(which(epi$refuge_scr < 0), 1))
p <- ggViolin(score="refuge_scr")
p1 <- p1 + xlab("") + ylab("pEPI")
p2 <- ggWrdBars(wrd_counts)
grid.arrange(p1, p2, ncol=2)

# TEST HERMIT
pstvs <- c('^sole', '^alone$', '^lone', '^single$', '^solitary$',
           '^unsociable', '^individual', '^hermit')
ngtvs <- c('^group', '^pack', '^altruistic', '^convivial',
           '^herd', '^flock', '^aggregate', '^swarm', '^horde', '^sociable')
wrd_counts <- senseScore(slt_nm="hermit_scr", pstvs, ngtvs)
if(!all(unique(epi$hermit_scr) %in% c(0, -1, 1, NA))) {
  stop('Unexpected result in sense score')
}
p1 <- ggViolin(score="hermit_scr")
p1 <- p1 + xlab("") + ylab("pEPI")
p2 <- ggWrdBars(wrd_counts)
grid.arrange(p1, p2, ncol=2)

# TEST SPECIALIST
pstvs <- c('^specialist$', '^unique$', '^specialised$', '^specialized$',
           '^idiosyncratic', '^adapted$', '^distinctive$', '^distinct',
           '^uncommon')
ngtvs <- c('^generalist', '^generic', '^unspecialised',
           '^unspecialized', '^opportunist$', '^opportunistic',
           '^omnivore$')
wrd_counts <- senseScore(slt_nm="spclst_scr", pstvs, ngtvs)
if(!all(unique(epi$spclst_scr) %in% c(0, -1, 1, NA))) {
  stop('Unexpected result in sense score')
}
p1 <- ggViolin(score="spclst_scr")
p1 <- p1 + xlab("") + ylab("pEPI")
p2 <- ggWrdBars(wrd_counts)
grid.arrange(p1, p2, ncol=2)

# TEST PRIMITIVE
pstvs <- c('^primitive', '^undeveloped', '^basic', '^simple', '^rudimentary$')
ngtvs <- c('^developed$', '^specialised$', '^specialized$', '^adapted$',
           '^advanced$', '^evolved$')
wrd_counts <- senseScore(slt_nm="prmtv_scr", pstvs, ngtvs)
if(!all(unique(epi$prmtv_scr) %in% c(0, -1, 1, NA))) {
  stop('Unexpected result in sense score')
}
p1 <- ggViolin(score="prmtv_scr")
p1 <- p1 + xlab("") + ylab("pEPI")
p2 <- ggWrdBars(wrd_counts)
grid.arrange(p1, p2, ncol=2)

# TEST ODDNESS
pstvs <- c('^odd', '^unusual', '^bizarre', '^aberrant',
           '^unconventional', '^peculiar', '^strange',
           '^atypical', '^unexpected', '^untypical',
           '^anomalous', '^remarkable')
ngtvs <- c('^standard', '^conventional')
wrd_counts <- senseScore(slt_nm="oddnss_scr", pstvs, ngtvs)
if(!all(unique(epi$oddnss_scr) %in% c(0, -1, 1, NA))) {
  stop('Unexpected result in sense score')
}
p1 <- ggViolin(score="oddnss_scr")
p1 <- p1 + xlab("") + ylab("pEPI")
p2 <- ggWrdBars(wrd_counts)
grid.arrange(p1, p2, ncol=2)

# CHECK MODELS
wilcox.test(epi$pepi[epi$refuge_scr == 1],
            epi$pepi[epi$refuge_scr == -1])
tapply(epi$pepi, factor(epi$refuge_scr), mean, na.rm=TRUE)
wilcox.test(epi$pepi[epi$hermit_scr == 1],
            epi$pepi[epi$hermit_scr == -1])
tapply(epi$pepi, factor(epi$hermit_scr), mean, na.rm=TRUE)
wilcox.test(epi$pepi[epi$spclst_scr == 1],
            epi$pepi[epi$spclst_scr == -1])
tapply(epi$pepi, factor(epi$spclst_scr), mean, na.rm=TRUE)
wilcox.test(epi$pepi[epi$prmtv_scr == 1],
            epi$pepi[epi$prmtv_scr == -1])
tapply(epi$pepi, factor(epi$prmtv_scr), mean, na.rm=TRUE)
wilcox.test(epi$pepi[epi$oddnss_scr == 1],
            epi$pepi[epi$oddnss_scr == -1])
tapply(epi$pepi, factor(epi$oddnss_scr), mean, na.rm=TRUE)

# ADDITIONAL
# random
all_wrds <- unlist(wrds)
all_wrds <- table(all_wrds)
rndm_wrds <- names(all_wrds > 200)  # find fairly common words
pstvs <- sample(rndm_wrds, 6)
ngtvs <- sample(rndm_wrds, 6)
wrd_counts <- senseScore(slt_nm="rndm_scr", pstvs, ngtvs)
wilcox.test(epi$pepi[epi$rndm_scr == 1],
            epi$pepi[epi$rndm_scr == -1])
tapply(epi$pepi, factor(epi$rndm_scr), mean, na.rm=TRUE)
# general
nwrds <- sapply(wrds, length)
plot(log(nwrds)~epi$pepi)  # no pattern of nwrds and pepi

# EPI TESTS
# refuge
# tapply(epi$epi, factor(epi$refuge_scr), mean, na.rm=TRUE)
# table(epi$refuge_scr)
# m1 <- aov(epi$epi ~ factor(epi$refuge_scr))
# m2 <- kruskal.test(epi$epi ~ factor(epi$refuge_scr))
# (m2)
# summary(m1)
# p <- ggViolin(score="refuge_scr", x="epi")
# p + xlab("Refuge concept") + ylab("EPI")

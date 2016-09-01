# INPUT
load(file.path('9_epi', 'res.RData'))
epi <- epi[!duplicated(epi$txid), ]
epi[['txid']] <- as.character(epi[['txid']])

epi$sstr_pepi <- NA
epi$sstr_cntrst_n <- NA
for(i in 1:nrow(epi)) {
  txid <- epi[['txid']][i]
  sstrs <- node_obj[[txid]][['sstr']]
  sstrs <- which(epi[['txid']] %in% sstrs)
  if(length(sstrs) > 0) {
    epi$sstr_pepi[i] <- mean(sapply(sstrs, function(x) epi[['pepi']][x]))
    epi$sstr_cntrst_n[i] <- mean(sapply(sstrs, function(x) epi[['cntrst_n']][x]))
  }
}

plot(log(epi$cntrst_n)~log(epi$sstr_cntrst_n))


write.csv(epi, file='~/Coding/Project-LF-Correlates/0_data/epi_scores.csv', row.names=FALSE)

plot(epi$pepi~epi$sstr_pepi)
which(epi[['nm']] == 'Monotremata')
which(epi[['nm']] == 'Theria')


cor.test(epi$pepi, epi$sstr_pepi, method='kendall')
cor.test(epi$pepi, epi$sstr_pepi, method='spearman')

pull <- epi$pepi < -10 & !is.na(epi$pepi)
epi[pull, 'nm']
abline(lm(epi$pepi ~ epi$sstr_pepi))
m1 <- lm(epi$pepi ~ epi$sstr_pepi)
epi <- epi[!is.na(epi$pepi), ]
epi$resids <- residuals(m1)
epi[['nm']][order(epi$resids)[1:100]]

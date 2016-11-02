### inc
library(plyr)
library(dplyr)

### par
batch_size <- 100
cores <- 60


### data
phen <- athaliana_phen(traits = "FRI", rows_order = "snp")
phen <- mutate(phen,
  log_FRI = log(FRI))

relmat <- athaliana_relmat()

gdat <- athaliana_snp()

### run GWAS
t1 <- proc.time()
gwas_lmer_wald <- assocLmer(log_FRI ~ (1|id), phen, relmat = list(id = relmat), data_snpcov = gdat, batch_size = batch_size, cores = cores, method = "Wald")
t2 <- proc.time()
(time_gwas_lmer_wald <- t2 - t1)

t1 <- proc.time()
gwas_lmer_lr <- assocLmer(FRI ~ (1|id), phen, relmat = list(id = relmat), data_snpcov = gdat, batch_size = batch_size, cores = cores, method = "LR")
t2 <- proc.time()
(time_gwas_lmer_lr <- t2 - t1)


### save
save(gwas_lmer_wald, time_gwas_lmer_wald, file = "gwas_lmer_wald.RData")
save(gwas_lmer_lr, time_gwas_lmer_lr, file = "gwas_lmer_lr.RData")

### inc
library(dplyr)
library(plyr)

library(lme4qtl)

### data
phen <- athaliana_phen(traits = "FRI")
relmat <- athaliana_relmat()

if(!exists("gdat")) {
  gdat <- athaliana_snp()
  
  gdat <- bind_cols(phen["FRI"], gdat)
}

dat <- gdat[1:20]

snps <- names(dat) %>% grep("^snp", ., value = TRUE)

mod0 <- relmatLmer(FRI ~ (1|id), dat, relmat = list(id = relmat), REML = FALSE)

#assoc <- assocLmer(FRI ~ (1|id), dat, relmat = list(id = relmat), data_snp_cov = subset(dat, select = snps), batch_size = 10)
stop()

### polygenic model
mod0 <- relmatLmer(FRI ~ (1|id), dat, relmat = list(id = relmat), REML = FALSE)

### testing SNPs for association
pvals <- laply(snps, function(x) {
  mod <- try({
    update(mod0, paste(". ~ . +", x))
  })

  if(class(mod)[1] != "try-error") {
    stats <- anova(mod, mod0)
    tab <- as.data.frame(stats)
    pval <- tab["mod", "Pr(>Chisq)"]
    stopifnot(!is.na(pval))
  }
  
  ifelse(class(mod)[1] == "try-error", as.numeric(NA), pval)
})


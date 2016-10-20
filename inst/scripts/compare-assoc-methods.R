### inc
library(rrBLUP)

library(dplyr)

### par
cores <- 10

### data
phen <- athaliana_phen(traits = "FRI", rows_order = "snp")
relmat <- athaliana_relmat()

gdat <- athaliana_snp(chr = 4)
gdat <- gdat[1:100]

### Method 1: rrBLUP::GWAS
# pheno: 
#  - Data frame where the first column is the line name (gid).
#  - Any column not designated as a fixed effect is assumed to be a phenotype.
#
# geno: 
#  - Data frame with the marker names in the first column.  
#  - The second and third columns contain the chromosome and map position
#  - Columns 4 and higher contain the marker scores for each line,
#  - The column names must match the line names in the "pheno" data frame.

pheno <- subset(phen, select = c("id", "FRI"))

ids <- gdat[[1]]
markers <- names(gdat)[-1]

gmat <- t(gdat[-1])
colnames(gmat) <- ids

geno <- bind_cols(data_frame(snp = markers, chr = 0, pos = 0), as_data_frame(gmat))

pheno <- as.data.frame(pheno)
geno <- as.data.frame(geno)
   
### run GWAS
assoc1_mm <- GWAS(pheno, geno, P3D = FALSE, n.core = cores, K = relmat, plot = FALSE)
assoc1_wm <- GWAS(pheno, geno, P3D = TRUE, n.core = cores, K = relmat, plot = FALSE)

### Method 2: lme4qtl
assoc2_wald <- assocLmer(FRI ~ (1|id), phen, relmat = list(id = relmat), data_snpcov = gdat, batch_size = 10, cores = cores, method = "Wald")
assoc2_lr <- assocLmer(FRI ~ (1|id), phen, relmat = list(id = relmat), data_snpcov = gdat, batch_size = 10, cores = cores, method = "LR")

### Ccompare results
tab <- left_join(
  data_frame(snp = assoc1_mm$snp, pval_mm = 10^(-assoc1_mm$FRI)),
  data_frame(snp = assoc1_wm$snp, pval_wm = 10^(-assoc1_wm$FRI)),
  by = "snp")

tab <- left_join(tab,
  data_frame(snp = assoc2_wald$snp, pval_wald = assoc2_wald$pval),
  by = "snp")
  
tab <- left_join(tab,
  data_frame(snp = assoc2_lr$snp, pval_lr = assoc2_lr$pval),
  by = "snp")  
  
  
  
  

#          snp    pval_mm    pval_wm    pval_wald      pval_lr
#        <chr>      <dbl>      <dbl>        <dbl>        <dbl>
#1  snp_124928 0.17972137 0.19215273 0.0090486009 0.0085826222

### inc
library(rrBLUP)

library(dplyr)

### data
phen <- athaliana_phen(traits = "FRI", rows_order = "snp")
relmat <- athaliana_relmat()

if(!exists("gdat")) {
  gdat <- athaliana_snp(chr = 4)
  gdat <- gdat[1:2]
}

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
#assoc1_mm <- GWAS(pheno, geno, P3D = FALSE, n.core = cores, K = relmat, plot = FALSE)
#assoc1_wm <- GWAS(pheno, geno, P3D = TRUE, n.core = cores, K = relmat, plot = FALSE)

### run GWAS by hand
ids <- pheno[["id"]]
y <- pheno[["FRI"]]
X <- matrix(1, length(y), 1)
g <- gdat[[2]]
K <- relmat

ind <- which(!is.na(y))
ids <- ids[ind]
y <- y[ind]
X <- X[ind, ]
g <- g[ind]
K <- K[ind, ind]

mod0 <- mixed.solve(y, K = K, X = X, return.Hinv = T)

y2 <- y
X3 <- cbind(X, g)

p <- ncol(X3)
v1 <- 1
v2 <- length(y2) - p

mod <- mixed.solve(y = y2, X = X3, K = K, return.Hinv = TRUE)
H2inv <- mod$Hinv

W <- crossprod(X3, H2inv %*% X3)
Winv <- solve(W)
beta <- Winv %*% crossprod(X3, H2inv %*% y2)
resid <- y2 - X3 %*% beta
s2 <- as.double(crossprod(resid, H2inv %*% resid))/v2
CovBeta <- s2 * Winv
Fstat <- beta[p]^2/CovBeta[p, p]
x <- v2/(v2 + v1 * Fstat)
pval <- pbeta(x, v2/2, v1/2)

### Raw data: `y` and `g`
assoc <- assocLmer(FRI ~ (1|id), phen, relmat = list(id = relmat), data_snpcov = gdat[1:2], method = "Wald")

# `lme4` model
dat <- data_frame(y = y, g = g, id = ids)
m <- relmatLmer(y ~ g + (1|id), dat, relmat = list(id = relmat))

# `lmmlite` model
library(lmmlite)

e <- eigen_rotation(K, y2, X3)
out <- fitLMM(e$Kva, e$y, e$X)

### Scaled data: `y` and `g`
# `lme4` model
dat <- data_frame(y = y, g = g, id = ids)
dat_sc <- mutate(dat, 
  y = scale(y), 
  g = scale(g))
m_sc <- relmatLmer(y ~ g + (1|id), dat_sc, relmat = list(id = relmat))

# `lmmlite` model
y2_sc <- scale(y2)
X3_sc <- X3
X3_sc[, 2] <- scale(X3_sc[, 2])

e_sc <- eigen_rotation(K, y2_sc, X3_sc)
out_sc <- fitLMM(e_sc$Kva, e_sc$y, e_sc$X)

### Compare two models: `m` and `m0`
dat <- data_frame(y = y, g = g, id = ids)

m0 <- relmatLmer(y ~ (1|id), dat, relmat = list(id = relmat))
m <- relmatLmer(y ~ g + (1|id), dat, relmat = list(id = relmat))

h2_m0 <- with(as.data.frame(VarCorr(m0)), vcov[1] / sum(vcov))
h2_m <- with(as.data.frame(VarCorr(m)), vcov[1] / sum(vcov))








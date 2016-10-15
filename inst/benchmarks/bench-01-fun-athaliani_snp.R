### inc
library(microbenchmark)

library(plyr)
library(ggplot2)

### par 
times <- 5

### small `n`
n <- 500

out1 <- microbenchmark(
  dplyr = athaliana_snp(method = "dplyr", n = n),
  fread = athaliana_snp(method = "fread", n = n),
  bycol = athaliana_snp(method = "bycol", n = n),
  times = times)
  
### moderate `n`
n <- 5e3

out2 <- microbenchmark(
  dplyr = athaliana_snp(method = "dplyr", n = n),
  bycol = athaliana_snp(method = "bycol", n = n),
  times = times)  
  
### relatevely large `n`
n <- 10e3

out3 <- microbenchmark(
  dplyr = athaliana_snp(method = "dplyr", n = n),
  bycol = athaliana_snp(method = "bycol", n = n),
  times = times)    

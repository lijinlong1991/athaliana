### inc
library(microbenchmark)

library(plyr)
library(ggplot2)

### par 
times <- 5

### small `n`
n <- 500

out <- microbenchmark(
  dplyr = athaliana_snp(method = "dplyr", n = n),
  fread = athaliana_snp(method = "fread", n = n),
  bycol = athaliana_snp(method = "bycol", n = n),
  times = times)
  
### moderate `n`
n <- 2e3

out <- microbenchmark(
  dplyr = athaliana_snp(method = "dplyr", n = n),
  bycol = athaliana_snp(method = "bycol", n = n),
  times = times)  

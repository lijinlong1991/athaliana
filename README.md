athaliana
=========

About
-----

`athaliana` is an R data package for the A. thaliana data set <https://github.com/Gregor-Mendel-Institute/atpolydb>.

Code examples
-------------

### Get necessary data with a few functions

Here we are inerested to run a polygenic model for one of the traits (`FRI`).

-   The necessary data are the table of phenotypes & the relatedness matrix (pre-calculated from the SNPs data, see `athaliana_compute_relmat` function).
-   The mixed model comes from the `lme4qtl` R package.

``` r
phen <- athaliana_phen(traits = "FRI")
relmat <- athaliana_relmat()

library(lme4qtl)
m <- relmatLmer(FRI ~ (1|id), phen, relmat = list(id = relmat))

m
```

    Linear mixed model fit by REML ['lmerMod']
    Formula: FRI ~ (1 | id)
       Data: phen
    REML criterion at convergence: 327.9857
    Random effects:
     Groups   Name        Std.Dev.
     id       (Intercept) 0.3232  
     Residual             0.4873  
    Number of obs: 164, groups:  id, 164
    Fixed Effects:
    (Intercept)  
          1.182

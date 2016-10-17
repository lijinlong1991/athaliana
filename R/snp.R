
#--------------------
# Relatedness matrix
#--------------------

#' @export
athaliana_relmat <- function(snp) 
{
  ### inc
  stopifnot(requireNamespace("rrBLUP"))
 
  A <- rrBLUP::A.mat(snp[-1])
}



#-----------------------------------
# Genetic Relatedness Matrix (GRM)
#-----------------------------------

#' @export
athaliana_compute_relmat <- function(snp) 
{
  ### prepare the matrix of genotypes: to be centered / scaled
  mat <- as.matrix(snp[-1])
  mat <- scale(mat, center = TRUE, scale = TRUE)

  ### var
  M <- ncol(mat)
  ids <- snp[["id"]]
  
  ### compute the var-covar matrix
  relmat <- tcrossprod(mat) / M
  
  rownames(relmat) <- ids
  colnames(relmat) <- ids  
  
  return(relmat)
}

#' @export
athaliana_rdata_relmat <- function() 
{
  "relmat.RData"
}

#' @export
athaliana_write_relmat <- function(relmat, 
  dir = file.path(athaliana_path(), athaliana_dir_rawdata()))
{
  file <- file.path(dir, athaliana_rdata_relmat())
  
  save(relmat, file = file)
}

#' @export
athaliana_relmat <- function(dir = file.path(athaliana_path(), athaliana_dir_rawdata()))
{
  file <- file.path(dir, athaliana_rdata_relmat())
  
  load(file)
  
  return(relmat)
}



#--------------------
# Files & Directories
#--------------------

#' @export
athaliana.guess.phendir <- function() 
{
  "~/git/variani/athaliana/rawdata/"
}

#' @export
athaliana.guess.phenfilename <- function() 
{
  url <- athaliana.guess.phenurl()
  basename(url)
}


#' @export
athaliana.guess.phenfile <- function() 
{
  file.path(athaliana.guess.phendir(), athaliana.guess.phenfilename())
}



#' @export
athaliana.guess.phenfile <- function() 
{
  file.path(athaliana.guess.phendir(), "phenotype_published_raw.tsv")
}

#--------------------
# URLs
#--------------------

#' @export
athaliana.guess.phenurl <- function() 
{
  "https://raw.githubusercontent.com/Gregor-Mendel-Institute/atpolydb/master/miscellaneous_data/phenotype_published_raw.tsv"
}


#--------------------
# Download
#--------------------

#' @export
athaliana.download.phen <- function(dir = athaliana.guess.phendir(), 
  url = athaliana.guess.phenurl(),
  ...)
{
  ### inc
  stopifnot(requireNamespace("utils"))
  
  ### args
  stopifnot(file.exists(dir))
  
  ### download 
  destfile <- file.path(dir, basename(url))
  ret <- download.file(url, destfile, ...)
}



#--------------------
# Files & Directories
#--------------------

#' @export
athaliana_dir_rawdata <- function() 
{
  "~/git/variani/athaliana/rawdata/"
}

#' @export
athaliana_filename_phen <- function() 
{
  url <- athaliana_url_phen()
  basename(url)
}


#' @export
athaliana_file_phen <- function() 
{
  file.path(athaliana_dir_rawdata(), athaliana_filename_phen())
}



#' @export
athaliana_file_phen <- function() 
{
  file.path(athaliana_dir_rawdata(), "phenotype_published_raw.tsv")
}

#--------------------
# URLs
#--------------------

#' @export
athaliana_url_phen <- function() 
{
  "https://raw.githubusercontent.com/Gregor-Mendel-Institute/atpolydb/master/miscellaneous_data/phenotype_published_raw.tsv"
}


#--------------------
# Download
#--------------------

#' @export
athaliana_download_phen <- function(dir = athaliana_dir_rawdata(), 
  url = athaliana_url_phen(),
  ...)
{
  ### inc
  stopifnot(requireNamespace("utils"))
  
  ### args
  if(!file.exists(dir)) {
    message(" * creating rawdata directory to store phenotypes: ", dir)
    ret <- dir.create(dir)
  }
  
  ### download 
  destfile <- file.path(dir, basename(url))
  ret <- utils::download.file(url, destfile, ...)
}


#----------------------------
# Load
#----------------------------

athaliana_phen <- function(dir = athaliana_dir_rawdata(),
  filename = athaliana_filename_phen(),
  traits, group, 
  names = c("clean", "raw"))
{
  ### arg
  names <- match.arg(names)
  
  ### inc
  stopifnot(requireNamespace("readr"))
  
  ### files/dirs
  stopifnot(file.exists(dir))
  
  file <- file.path(dir, filename)
  stopifnot(file.exists(file))
  
  ### read
  phen <- readr::read_tsv(file)
  
  ### names
  if(names == "clean") {
    nms <- names(phen)
    
    nms <- str_replace(nms, "^[0-9]*_", "")

    nms <- str_replace_all(nms, " ", "_")
    nms <- str_replace_all(nms, "-", "_")
    
    nms <- str_replace(nms, "<i>", "")
    nms <- str_replace(nms, "</i>", "")
    
    names(phen) <- nms
  }
  
  ### select columns  
  if(!all(c(missing(traits), missing(group)))) {
    vars <- phen %>% names %>% head(2)
    
    if(!missing(group)) {
      traits <- switch(group,
        "strong" = athaliana_traits_strong(),
        stop())
      stopifnot(all(traits %in% names(phen)))
      
      vars <- c(vars, traits)
    }
    
    if(!missing(traits)) {
      stopifnot(all(traits %in% names(phen)))
      
      vars <- c(vars, traits)
      vars <- unique(vars)
    }
    
    phen <- select(phen, one_of(vars))
  }
    
  return(phen)
}



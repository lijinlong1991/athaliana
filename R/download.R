
#--------------------
# Files & Directories
#--------------------

#' @export
athaliana_path <- function() 
{
  "~/git/variani/athaliana/" %>% path.expand
}


#' @export
athaliana_dir_rawdata <- function() 
{
  "rawdata"
}

#----------------------------
# Phen Files & Directories
#----------------------------

#' @export
athaliana_filename_phen <- function() 
{
  url <- athaliana_url_phen()
  basename(url)
}

#' @export
athaliana_file_phen <- function() 
{
  file.path(athaliana_path(), athaliana_dir_rawdata(), athaliana_filename_phen())
}

#----------------------------------
# Genotypes Files & Directories
#----------------------------------

#' @export
athaliana_filename_snp <- function() 
{
  "call_method_32.b"
}

#' @export
athaliana_path_snp <- function() 
{
  file.path(athaliana_path(), 
    athaliana_dir_rawdata(),
    "Network", "Data", "250k", "db", "dataset")
}


#' @export
athaliana_file_snp <- function() 
{
  file.path(athaliana_path_snp(), 
    athaliana_filename_snp())
}

#--------------------
# URLs
#--------------------

#' @export
athaliana_url_phen <- function() 
{
  "https://raw.githubusercontent.com/Gregor-Mendel-Institute/atpolydb/master/miscellaneous_data/phenotype_published_raw.tsv"
}

#' @export
athaliana_url_snp <- function() 
{
  "https://github.com/Gregor-Mendel-Institute/atpolydb/raw/master/250k_snp_data/call_method_32.tar.gz"
}

#--------------------
# Download
#--------------------

#' @export
athaliana_download_phen <- function(dir = file.path(athaliana_path(), athaliana_dir_rawdata()), 
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

#' @export
athaliana_download_snp <- function(dir = file.path(athaliana_path(), athaliana_dir_rawdata()), 
  url = athaliana_url_snp(),
  extract = TRUE,
  ...)
{
  ### inc
  stopifnot(requireNamespace("utils"))
  
  ### args
  if(!file.exists(dir)) {
    message(" * creating rawdata directory to store snps: ", dir)
    ret <- dir.create(dir)
  }
  
  ### download 
  destfile <- file.path(dir, basename(url))
  ret <- utils::download.file(url, destfile, ...)

  ### extract  
  if(extract) {
    untar(destfile, exdir = dir)
  }
}


#----------------------------
# Load
#----------------------------

athaliana_phen <- function(file = file.path(athaliana_path(), athaliana_dir_rawdata(), athaliana_filename_phen()),
  traits, group, 
  names = c("clean", "raw"))
{
  ### arg
  names <- match.arg(names)
  
  ### inc
  stopifnot(requireNamespace("readr"))
  
  ### files/dirs
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
    
    phen <- subset(phen, select = vars)
  }
    
  return(phen)
}


#' @export
athaliana_snp <- function(file = athaliana_file_snp())
{
  ### inc
  stopifnot(requireNamespace("readr"))
  
  ### files/dirs
  stopifnot(file.exists(file))
  
  ### read: version 1 (slow)
  #dat <- readr::read_csv(file, skip = 1)

  #annot <- select(dat, 1:2)
  #dat <- t(select(dat, -c(1:2))) 
  
  ### read: version 2 via reading line by line
  # @ http://stackoverflow.com/questions/17288197/reading-a-csv-file-organized-horizontally

  # read_lines(file, skip = 3, 1) %>% str_split(",") %>% unlist
}

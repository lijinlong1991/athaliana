
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

#' @export
athaliana_feather_snp <- function() 
{
  "snp.feather"
}

#' @export
athaliana_feather_annot <- function() 
{
  "annot.feather"
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


#--------------------------------------
# Read/feather SNP annotation data
#--------------------------------------

#' @export
athaliana_write_annot <- function(dir = file.path(athaliana_path(), athaliana_dir_rawdata()), 
  ...)
{
  ### inc
  stopifnot(requireNamespace("feather"))

  ### read SNP annot. data
  annot <- athaliana_read_annot(...)

  ### read feather
  path <- file.path(dir, athaliana_feather_annot())
  feather::write_feather(annot, path) 
} 


#' @export
athaliana_read_annot <- function(file = athaliana_file_snp(), 
  verbose = 1)
{
  annot <- read_csv(file, skip = 1, 
    col_type = cols_only(Chromosome = col_integer(), Positions = col_integer()))
  
  annot <- bind_cols(tibble(snp = paste0("snp_", 1:nrow(annot))), annot)
  
  names(annot) <- c("snp", "chr", "pos")
  
  return(annot)
}

#' @export
athaliana_annot <- function(dir = file.path(athaliana_path(), athaliana_dir_rawdata()), 
  ...)
{
  ### inc
  stopifnot(requireNamespace("feather"))

  ### write feather
  path <- file.path(dir, athaliana_feather_annot())
  feather::read_feather(path)
} 

#--------------------------------------
# Read/feather SNP data
#--------------------------------------

#' @export
athaliana_write_snp <- function(dir = file.path(athaliana_path(), athaliana_dir_rawdata()), 
  ...)
{
  ### inc
  stopifnot(requireNamespace("feather"))

  ### read SNP data
  snp <- athaliana_read_snp(...)

  ### read feather
  path <- file.path(dir, athaliana_feather_snp())
  feather::write_feather(snp, path) 
} 

#' @export
athaliana_snp <- function(dir = file.path(athaliana_path(), athaliana_dir_rawdata()), 
  chr,
  ...)
{
  ### inc
  stopifnot(requireNamespace("feather"))

  ### write feather
  path <- file.path(dir, athaliana_feather_snp())
  snp <- feather::read_feather(path)
  
  ### filter by chr
  if(!missing(chr)) {
    annot <- athaliana_annot(dir = dir)
    
    chr_val <- chr
    snps <- with(annot, snp[chr %in% chr_val])
    
    snp <- subset(snp, select = c("id", snps))
  }
  
  ### return
  return(snp)
} 


#' @export
athaliana_read_snp <- function(file = athaliana_file_snp(), 
  n, method = c("dplyr", "fread", "bycol"),
  format = c("numeric", "raw"),
  verbose = 1)
{
  ### args
  method <- match.arg(method)
  format <- match.arg(format)  
  
  missing.n <- missing(n)
  
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
  # code example: read_lines(file, skip = 3, 1) %>% str_split(",") %>% unlist
  
  # head line
  ids <- athaliana_ids_snp(file = file)
  
  f <- function(x) 2 * (x - 1.5)
  
  ### read genotypes by means of `method`
  if(method == "dplyr") { 
    if(missing.n) {
      lines <- readr::read_lines(file, skip = 2)
    } else {
      lines <- readr::read_lines(file, skip = 2, n_max = n)
    }
    
    mat <- switch(format, 
      "numeric" = {
        lines %>%
        lapply(. %>% str_split(",") %>% unlist %>% tail(-2) %>%
          as.factor %>% as.numeric %>% f %>%
          tibble) %>%
        do.call("bind_cols", .)
      },
      "raw" = {
        lines %>%
        lapply(. %>% str_split(",") %>% unlist %>% tail(-2) %>%
          tibble) %>%
        do.call("bind_cols", .)
      },
      stop())
    
    # add names of SNPs   
    names(mat) <- paste0("snp_", 1:ncol(mat))
    
    # add columns of IDs
    mat <- bind_cols(tibble(id = ids), mat)
  } 
  else if(method == "bycol") {
    lines <- readr::read_lines(file, skip = 2, n_max = n)
    
    mat <- matrix(as.numeric(NA), nrow = length(ids), ncol = length(lines))
    for(i in 1:length(lines)) {
      if(verbose) {
        if(i %% 1e4 == 0) {
          cat(" * marker", i, "/", length(lines), "\n")
        }
      }
  
      val <- lines[i] %>% str_split(",") %>% unlist %>% tail(-2) %>%
        as.factor %>% as.numeric %>% f
    
      mat[, i] <- val
    }
    rownames(mat) <- ids  
  } 
  else if(method == "fread") {
    requireNamespace("data.table")
    
    dat <- data.table::fread(file, skip = 1, header = TRUE, nrows = n)
    
    mat <- matrix(as.numeric(NA), nrow = length(ids), ncol = nrow(dat))
    for(i in 1:nrow(dat)) {
      if(verbose) {
        if(i %% 1e4 == 0) {
          cat(" * marker", i, "/", nrow(dat), "\n")
        }
      }
      
      val <- as.character(dat[i, ]) %>% tail(-2) %>%
        as.factor %>% as.numeric %>% f
      
      mat[, i] <- val
    }
  } 
  else {
    stop("`method` unknown")
  }

  return(mat) 
}

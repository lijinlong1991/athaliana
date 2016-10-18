#---------------------------------------
# Get phen from a (downloaded) file 
#---------------------------------------

athaliana_phen <- function(file = file.path(athaliana_path(), athaliana_dir_rawdata(), athaliana_filename_phen()),
  traits, group, 
  names = c("clean", "raw"),
  rows_order = c("none", "snp"))
{
  ### arg
  names <- match.arg(names)
  rows_order <- match.arg(rows_order)
    
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
  
  ### derive new columns
  phen <- mutate(phen, 
    id = as.character(ecotype_id))

  ### sort rows 
  if(rows_order == "snp") {
    # @ http://stackoverflow.com/questions/26548495/reorder-rows-using-custom-order
    ids_snp <- athaliana_ids_snp()
    ids_phen <- phen$id
    
    stopifnot(length(ids_snp) == length(ids_phen))
    stopifnot(all(ids_snp %in% ids_phen))
    stopifnot(all(ids_phen %in% ids_snp))    
    
    ids_fac <- factor(ids_phen, levels = ids_snp)
    ids_ind <- order(ids_fac)
    
    # arrange `phen`
    phen <- phen[ids_ind, ]
    
    stopifnot(all(phen$id == ids_snp))
  }
  
  return(phen)
}

#----------------------------
# IDs
#----------------------------

#' @export
athaliana_ids_phen <- function(...)
{
  phen <- athaliana_phen(..., traits = NULL)

  phen[["id"]]
}

#' @export
athaliana_ids_snp <- function(file = athaliana_file_snp())
{
  ### inc
  stopifnot(requireNamespace("readr"))
  
  ### files/dirs
  stopifnot(file.exists(file))
  
  # head line
  hline <- readr::read_lines(file, skip = 1, n_max = 1)
  ids <- hline %>% str_split(",") %>% unlist %>% tail(-2)
  
  return(ids)
}

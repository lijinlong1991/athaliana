
#--------------------
# NA's
#--------------------

#' @export
athaliana_stat_phen_na <- function(phen)
{
  ### args
  stopifnot(!missing(phen))
  
  ### compute stats
  nobs <- nrow(phen)

  num_na <- phen %>% summarize_each(funs(. %>% is.na %>% sum)) %>% as.integer

  data_frame(variable = names(phen), num_na = num_na) %>%
    arrange(desc(num_na)) %>%
    mutate( 
      variable = factor(variable, levels = unique(variable)),
      num_obs = nobs - num_na)
}

#' @export
athaliana_plot_phen_na <- function(phen)
{
  ### args
  stopifnot(!missing(phen))
  
  ### compute stats
  nobs <- nrow(phen)
  
  tab <- athaliana_stat_phen_na(phen)
  
  ggplot(tab, aes(variable, num_obs)) + 
    geom_point() + 
    ylim(c(0, nobs)) + 
    geom_hline(yintercept = c(nobs, nobs/2), linetype = c(3, 3)) + 
    coord_flip() + 
    theme_linedraw()
}

#' Traits with strong peaks of association
#'
#' \code{FT_LD}: flowering time under long days;
#' \code{Dormancy} / code{8_}: seed dormancy;
#' \code{avrRpm}: disease resistance;
#' \code{FRI}: FRI gene expression.
#'
#' @references \url{http://potatobreeding.cals.wisc.edu/wp-content/uploads/sites/21/2014/01/GWAS_tutorial.pdf}
#'
#' @export
athaliana_traits_strong <- function()
{
  c("LD", "Seed_Dormancy", "avrRpm1", "FRI")
}

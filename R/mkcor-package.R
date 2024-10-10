#' mkcor
#'
#' @description A package for fitting the correlated Mk model to character data on a phylogeny
#'
#' @section Functions for model fitting and simulation:
#' \itemize{
#'   \item \code{\link{mkcor_em}}
#'   \item \code{\link{mkcor_loglik}}
#'   \item \code{\link{mkcor_fit_em}}
#'   \item \code{\link{mkcor_fit}}
#'   \item \code{\link{simulate_mkcor}}
#' }
#' @section Functions for working with trees:
#' \itemize{
#'   \item \code{\link{as.ephylo}}
#'   \item \code{\link{ephylo}}
#'   \item \code{\link{ephylo_preorder}}
#'   \item \code{\link{ephylo_postorder}}
#'   \item \code{\link{ephylo_tips}}
#'   \item \code{\link{ephylo_rclade}}
#'   \item \code{\link{read_newick}}
#' }
#' @section Vignettes:
#' \code{vignette("CorrelatedMk")}
#' @docType package
#' @name mkcor
"_PACKAGE"
#' Metabolite names and IDs.
#'
#' A dataset containing metabolite names and IDs from Recon 1 and 2.2.
#' These are also stored in the correct order to reproduce Fig 3 from the Recon 2 article.
#' 
#' @format A data frame with 54 rows and 3 variables:
#' \describe{
#'     \item{Metab.names}{The metabolite names}
#'     \item{Metab.Recon1.IDs}{Metabolite IDs from Recon 1}
#'     \item{Metab.Recon2v2.IDs}{Metabolite IDs from Recon 2.2}
#'     \item{Metab.EX.IDs}{Metabolite IDs in the EX_xxx(e) format}
#' }
"metab_names"

#' IEM names and IDs.
#'
#' A dataset containing IEM names and IDs from Fig3 of the Recon 2 article.
#' These are also stored in the correct order to reproduce Fig 3 from the Recon 2 article.
#' 
#' @format A data frame with 49 rows and 4 variables:
#' \describe{
#'     \item{IEM.full}{The full IEM names, may contain lower cases}
#'     \item{IEM.code}{IEM codes used in sampling for identification}
#'     \item{IEM.full.upper}{The full IEM names with consistent capitalisation}
#'     \item{IEM.short}{Short versions of the IEM names for figures}
#' }
"iem_names"

#' Observed metabolite changes in IEM patients.
#'
#' A dataset containing the observed metabolite changes in patients from Fig3 of the Recon 2 article.
#' 
#' @format A data frame with 54 rows and 50 variables:
#' \describe{
#'     \item{Metab}{All 54 metabolites from Fig 3 of the Recon 2 article}
#'     \item{IEM names}{All the other columns are the 49 IEMs from Fig 3 of the Recon 2 article}
#' }
"obs"
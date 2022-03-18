# Functions for calculating


#' Calculate the difference between random pairs of flux samples
#' 
#' This function uses sdata to calculate random subtractions between c2 and c1 sampling results.
#' 
#' @param sdata An sdata list containing c1 and c2 dataframes
#' @param seed A seed to set for reproducibility
#' @param divfactor A division factor to divide the number of samples by for less intensive calculations
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace_all
#' @importFrom dplyr select bind_cols rename_with mutate group_by summarise
#' @importFrom tidyr gather
#' @returns A "long" dataframe containing all diff samples for all metabolites
#' @export
calc_diff = function(sdata, seed = NA, divfactor=1) {
  if (!(is.na(seed))){
    set.seed = seed
  }
  n = dim(sdata$MUT)[1]/divfactor
  sdata$MUT = sdata$MUT %>% reshape2::melt(variable.name="Metab") %>% mutate("Type" = "MUT")
  sdata$WT = sdata$WT %>% reshape2::melt(variable.name="Metab") %>% mutate("Type" = "WT")
  gc()
# 
#   sdata_small = list(MUT = sdata$MUT[1:1000,], WT = sdata$WT[1:1000,])
#   n = dim(sdata_small$MUT)[1]/divfactor
#   sdata_small$MUT = sdata_small$MUT %>% reshape2::melt(variable.name="Metab") %>% mutate("Type" = "MUT")
#   sdata_small$WT = sdata_small$WT %>% reshape2::melt(variable.name="Metab") %>% mutate("Type" = "WT")

  sdata = data.table::rbindlist(sdata)
  gc()
  diff = sdata %>% group_by(Metab) %>%
    # Difference: mut - wt (c2 - c1)
    summarise(Value = sample(value[Type=="MUT"], n) - sample(value[Type=="WT"], n))

return(diff)
}

#' Calculate a pseudo-z-score of a difference distribution
#' 
#' This function uses a calculated diff (calc_diff) to calculate a pseudo-z-score for the 
#' distribution.
#' 
#' @param diff A "long" dataframe containing Metabs and diff values
#' @param disease.code The 3-letter disease code for the current disease
#' @returns A zscore dataframe containing zscores for each metabolite
#' @importFrom magrittr %>%
#' @export
calc_zscore = function(diff){
  zscore = diff %>%
    group_by(Metab) %>%
    summarise("zscore" = mean(Value) / sd(Value))
  return(zscore)
}
# Functions for calculating


#' Calculate the difference between random pairs of flux samples
#' 
#' This function uses sdata to calculate random subtractions between c2 and c1 sampling results.
#' 
#' @param sdata An sdata list containing WT and MUT dataframes
#' @param seed A seed to set for reproducibility
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace_all
#' @importFrom dplyr select bind_cols rename_with mutate group_by summarise
#' @importFrom tidyr gather
#' @import bigmemory
#' @returns A "long" dataframe containing all diff samples for all metabolites
#' @export
calc_diff = function(sdata, seed = NA) {
  if (!(is.na(seed))){
    set.seed = seed
  }
  #n = dim(sdata$MUT)[1]/divfactor
  #sdata$MUT = sdata$MUT %>% reshape2::melt(variable.name="Metab") %>% mutate("Type" = "MUT")
  #sdata$WT = sdata$WT %>% reshape2::melt(variable.name="Metab") %>% mutate("Type" = "WT")
  
  metab_names = names(sdata$WT)
  MUT_matrix = as.big.matrix(sdata$MUT)
  rand_rows_MUT = sample(nrow(MUT_matrix))
  WT_matrix = as.big.matrix(sdata$WT)
  rand_rows_WT = sample(nrow(WT_matrix))
  diff = as.data.frame(MUT_matrix[rand_rows_MUT, ] - WT_matrix[rand_rows_WT,])
  colnames(diff) = metab_names
  gc()
#   sdata_small = list(MUT = sdata$MUT[1:1000,], WT = sdata$WT[1:1000,])
#   n = dim(sdata_small$MUT)[1]/divfactor
#   sdata_small$MUT = sdata_small$MUT %>% reshape2::melt(variable.name="Metab") %>% mutate("Type" = "MUT")
#   sdata_small$WT = sdata_small$WT %>% reshape2::melt(variable.name="Metab") %>% mutate("Type" = "WT")

  # sdata = data.table::rbindlist(sdata)
  # gc()
  # diff = sdata %>% group_by(Metab) %>%
  #   # Difference: mut - wt (c2 - c1)
  #   summarise(Value = sample(value[Type=="MUT"], n) - sample(value[Type=="WT"], n))

return(diff)
}

#' Calculate a pseudo-z-score of a difference distribution
#' 
#' This function uses a calculated diff (calc_diff) to calculate a pseudo-z-score for the 
#' distribution.
#' 
#' @param diff A "long" dataframe containing Metabs and diff values
#' @returns A zscore dataframe containing zscores for each metabolite
#' @importFrom magrittr %>%
#' @export
calc_zscore = function(diff){
  zscore = apply(diff, 2, function(x) mean(x)/sd(x)) %>%
    as.data.frame() %>%
    rename("zscore" = ".") %>%
    rownames_to_column("Metab")
    
  # zscore = diff %>%
  #   group_by(Metab) %>%
  #   summarise("zscore" = mean(Value) / sd(Value))
  return(zscore)
}


#' Calculate the mean of each exchange reaction in sdata
#' 
#' This function calculates the mean for each column (=exchange reaction = metabolite) and transforms the data into a
#' format usable by plot_scatter.
#' 
#' @param sdata An sdata list containing WT and MUT dataframes
#' @returns A means list containing WT and MUT mean dataframes
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
calc_sdata_means = function(sdata){
  # Calculate averages for each metabolite's flux
  avg_WT = sdata$WT %>%
    summarise(across(everything(), mean))
  avg_MUT = sdata$MUT %>%
    summarise(across(everything(), mean))
  # Pivot the dataframes so rows are metabolites and means are in a column
  means_WT = avg_WT %>% t() %>% as.data.frame() %>%
    rownames_to_column("Metab") %>%
    rename("mean_WT" = V1)
  means_MUT = avg_MUT %>% t() %>% as.data.frame() %>%
    rownames_to_column("Metab") %>%
    rename("mean_MUT" = V1)
  means = list(means_WT, means_MUT)
  names(means) = c("WT","MUT")
  return(means)
}


#' Calculate the densities of each metabolite in sdata
#' 
#' This function calculates the density for each column (=exchange reaction = metabolite).
#' 
#' @param sdata An sdata list containing WT and MUT dataframes
#' @returns A means list containing WT and MUT mean dataframes
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
calc_density = function(sdata){
  # Calculate densities for each metabolite's flux
  densities = list(apply(sdata$WT, 2, density), apply(sdata$MUT, 2, density))
  names(densities) = c("WT","MUT")
  return(densities)
}


#' Filter the metabolites based on their zscore
#' 
#' This function filters metabolites based on their previously calculated zscore using sdata.
#' 
#' @param sdata An sdata list containing WT and MUT dataframes
#' @returns A means list containing WT and MUT mean dataframes
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
calc_filter_sdata = function(sdata, zscores, thr=2, max_n=10){
  metabs_filtered = zscores %>% arrange(desc(abs(zscore))) %>%
    slice_head(n = max_n) %>% filter(abs(zscore) > thr)
  sdata_all_filtered = lapply(sdata, function(x) x %>% select(metabs_filtered$Metab))
  return(sdata_all_filtered)
}


#' Filter the metabolites based on their zscore
#' 
#' This function filters metabolites based on their previously calculated zscore using densities.
#' 
#' @param sdata An sdata list containing WT and MUT dataframes
#' @returns A means list containing WT and MUT mean dataframes
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
calc_filter_density = function(sdata, zscores, thr=2, max_n=10){
  metabs_filtered = zscores %>% arrange(desc(abs(zscore))) %>%
    slice_head(n = max_n) %>% filter(abs(zscore) > thr)
  densities = calc_density(sdata)
  d_all_filtered = lapply(densities, function(x) x)
  d_all_filtered = list(densities$WT[names(densities$WT) %in% metabs_filtered$Metab],
                             densities$MUT[names(densities$MUT) %in% metabs_filtered$Metab])
  return(d_all_filtered)
}


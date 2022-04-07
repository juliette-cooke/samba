#' Filter the metabolites based on their zscore
#' 
#' This function filters metabolites based on their previously calculated zscore and density.
#' 
#' @param d_all A list containing WT and MUT dataframes
#' @returns A means list containing WT and MUT mean dataframes
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
filter_all_density = function(d_all, zscores, thr=2, max_n=10){
  d_all_filtered = list()
  for (i in 1:length(d_all)){
    metabs_filtered = zscores[[i]] %>% arrange(desc(abs(zscore))) %>%
      slice_head(n = max_n) %>% filter(abs(zscore) > thr)
    d_all_filtered[[i]] = list(d_all[[i]]$WT[names(d_all[[i]]$WT) %in% metabs_filtered$Metab],
                          d_all[[i]]$MUT[names(d_all[[i]]$MUT) %in% metabs_filtered$Metab])
    names(d_all_filtered[[i]]) = c("WT","MUT")
    names(d_all_filtered)[i] = names(d_all)[i]
  }
  return(d_all_filtered)
}
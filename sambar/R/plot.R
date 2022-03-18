# Functions for plotting

#' Plot distributions
#' 
#' Compare predicted and observed for any set of metabolites and KOs. You must provide both the pred and obs datasets.
#' The first column of both obs and pred must be the "Metab" column.
#' The other columns are the condition columns. The column names must match between obs and pred.
#' Obs must contain "+1", "-1" and "00". Pred must be numeric.
#' Conditions (column names) must match the ones in condition_order.
#' Metabolites (first column) must match the ones in metab_order.
#' 
#' @param pred The predicted dataframe with Metabs as rows and conditions as columns, with Metabs as column 1
#' @param obs The observed dataframe with Metabs as rows and conditions as columns, with Metabs as column 1
#' @param title The title to be used for the plot
#' @param metab_order The order in which to plot the metabolites (must match metabolites in Metabs)
#' @param condition_order The order in which to plot the onditions (must match column names)
#' @param squash Whether to squash extreme values (squash>0) or not (squash=0)
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join slice_head filter arrange select mutate 
#' @importFrom tidyr gather pivot_longer
#' @import ggplot2
#' @export
plot_distrib = function(sdata, zscore, thr, max_n=10, case="Disease"){
  thr = 2
  max_n = 10
  top_metab_zscore = zscore %>% arrange(desc(abs(zscore))) %>%
    slice_head(n = max_n) %>% filter(abs(zscore) > thr)
  
  sdata_filter = rbind(sdata$WT %>% select(top_metab_zscore$Metab) %>% mutate(Type = "WT"),
                       sdata$MUT %>% select(top_metab_zscore$Metab) %>% mutate(Type = "MUT")) %>%
    pivot_longer(-Type, names_to = "Metab")

  sdata_filter %>%
    ggplot(aes(x=value, fill=Type))+
    geom_density(alpha = 0.5, colour=NA)+
    facet_wrap(~Metab,nrow = 1, scales = "free")+
    theme_minimal()+
    ggtitle(paste0(paste0("Top ",max_n," flux sampling distributions filtered with a z-score threshold of ",thr)))+
    ylab(case)+
    xlab("Flux value")
}
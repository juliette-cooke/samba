# Functions for plotting

#' Plot distributions
#' 
#' Compare WT and MUT flux distributions for the top set of zscore metabolites.
#' 
#' @param sdata A list containing a WT dataframe and a MUT dataframe. Each dataframe contains X samples for metabolites
#' @param zscore A dataframe of zscore calculated using calc_zscore
#' @param thr The threshold for filtering zscores
#' @param max_n The max number of metabolites to show
#' @param case The name of the disease or condition
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join slice_head filter arrange select mutate 
#' @importFrom tidyr gather pivot_longer
#' @import ggplot2
#' @export
plot_distrib = function(sdata, zscore, thr, max_n=10, case="Disease"){
  #thr = 2
  #max_n = 10
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


#' Plot a scatter plot of WT means vs MUT means
#' 
#' Plots a scatter plot with WT flux means on the x-axis and MUT flux means on the y-axis. The points are coloured
#' by zscore.
#' 
#' @param sdata A list containing a WT dataframe and a MUT dataframe. Each dataframe contains X samples for metabolites
#' @param zscore A dataframe of zscore calculated using calc_zscore
#' @param thr The threshold for filtering zscores
#' @param case The name of the disease or condition
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom tidyr gather pivot_longer
#' @importFrom ggforce facet_zoom
#' @importFrom tibble rownames_to_column
#' @importFrom viridis scale_colour_viridis
#' @import ggplot2
#' @export
plot_scatter = function(sdata, zscore, thr, case="Disease"){
  # Calculate averages for each metabolite's flux
  avg_WT = sdata$WT %>%
    summarise(across(everything(), mean))
  avg_MUT = sdata$MUT %>%
    summarise(across(everything(), mean))
  # Pivot the dataframes so rows are metabolites and means are in a column
  t1 = avg_WT %>% t() %>% as.data.frame() %>%
    rownames_to_column("Metab") %>%
    rename("mean_WT" = V1)
  t2 = avg_MUT %>% t() %>% as.data.frame() %>%
    rownames_to_column("Metab") %>%
    rename("mean_MUT" = V1)
  # Join the zscore with the means
  z_m = zscore %>% left_join(., t1 , by="Metab") %>% left_join(., t2 , by="Metab")
  
  #thr = 1
  # Filter the zscore based on the provided threshold
  z_m_filtered = z_m %>% filter(zscore > thr | zscore < -thr)
  # Create the breaks for the legend
  z_breaks = round(seq(from=min(z_m_filtered$zscore), to=max(z_m_filtered$zscore),by=5))
  z_breaks = z_breaks[abs(z_breaks-thr) > 4]
  z_breaks_thr = c(z_breaks, thr, -thr)
  z_labels = c(z_breaks, paste0(" thr=",thr), paste0("-thr=", -thr))
  
  z_m_filtered %>% 
    ggplot(aes(x=mean_WT, y=mean_MUT, colour=zscore))+
    geom_point(size = 4, alpha = 0.8)+
    geom_abline(intercept = 0, linetype = "dashed")+
    # coord_cartesian(xlim = c(min(min(z_m_filtered$mean_WT),min(z_m_filtered$mean_MUT)), 
    #        max(max(z_m_filtered$mean_WT),max(z_m_filtered$mean_MUT))),
    #        ylim = c(min(min(z_m_filtered$mean_WT),min(z_m_filtered$mean_MUT)), 
    #        max(max(z_m_filtered$mean_WT),max(z_m_filtered$mean_MUT))))+
    facet_zoom(xlim = c(-1,1), ylim = c(-1,1))+
    scale_colour_viridis(option = "magma", breaks = z_breaks_thr, labels= z_labels)+
    theme_light()+
    ggtitle(paste0("WT mean fluxes vs KO mean fluxes for ", case))+
    xlab("WT mean flux")+
    ylab("KO mean flux")
  
}
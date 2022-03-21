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


#' Plot multiple distributions
#' 
#' Compare WT and MUT flux distributions for the top set of zscore metabolites for multiple diseases.
#' 
#' @param sdata A list containing a WT dataframe and a MUT dataframe. Each dataframe contains X samples for metabolites
#' @param zscore A dataframe of zscore calculated using calc_zscore
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate 
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @importFrom ggh4x facet_nested_wrap
#' @export
plot_multi_distrib = function(sdata_filtered, thr, max_n=10){
  #thr = 2
  #max_n = 10
  #sdata_filtered = sdata_all_filtered
  i=1
  sdata_filter = data.frame()
  n = length(sdata_filtered)
  for(i in 1:n){
    temp = rbind(sdata_filtered[[i]]$WT %>% mutate(Type = "WT"),
                 sdata_filtered[[i]]$MUT %>% mutate(Type = "MUT")) %>%
      pivot_longer(-Type, names_to = "Metab") %>%
      mutate(Case = names(sdata_filtered)[i])
    sdata_filter = rbind(sdata_filter, temp)
  }
  
  sdata_filter %>%
    ggplot(aes(x=value, fill=Type))+
    geom_density(alpha = 0.5, colour=NA)+
    facet_nested_wrap(~Case + Metab, nrow = n, scales = "free")+
    theme_minimal()+
    ggtitle(paste0(paste0("Top flux sampling distributions filtered with a z-score threshold of ",thr)))+
    ylab("")+
    xlab("Flux value")
}


#' Plot a basic demo distribution
#' 
#' Compare WT and MUT flux distributions for the top set of zscore metabolites for multiple diseases.
#' 
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export
plot_basic_distrib = function(){
  example_d = data.frame(WT=rnorm(5000, 0, 1), MUT=rnorm(5000, 4, 1))
  example_d %>% pivot_longer(cols = c(WT,MUT),names_to = "Type", values_to = "value") %>%
    ggplot(aes(x=value, fill=Type))+
    geom_density(alpha = 0.5, colour=NA)+
    theme_minimal()+
    theme(plot.margin = unit(c(0,0,0,0), "cm"))+
    ggtitle(paste0("Example flux sampling distributions for metabolite X"))+
    ylab("")+
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
plot_scatter = function(sdata, zscore, means,thr, case="Disease"){
  # Join the zscore with the means
  z_m = zscore %>% left_join(., means$WT , by="Metab") %>% left_join(., means$MUT , by="Metab")
  
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
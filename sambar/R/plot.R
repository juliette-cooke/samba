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
    # facet_wrap(~Metab,nrow = 1, scales = "free")+
    facet_wrap(~Metab,nrow = 1, scales = "free_y")+
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
  #i=1
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


#' Plot multiple density distributions
#' 
#' Compare WT and MUT flux distributions for the top set of zscore metabolites for multiple diseases using filtered 
#' densities as inputs instead of raw sdata.
#' 
#' @param densities_filt A list containing a WT list and a MUT list. Each list contains X densities for metabolites
#' @param zscore A dataframe of zscore calculated using calc_zscore
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @importFrom ggh4x facet_nested_wrap
#' @export
plot_multi_density_distrib = function(densities_filt, metab_dict, thr=2, max_n=10){
  #thr = 2
  #max_n = 10
  #densities_filt = d_all_filtered
  #i=1
  d_filter = list()
  n = length(densities_filt)
  for(i in 1:n){
    # Extract x and y from the density objects and pivot into long format 
    temp = rbind(cbind(as.data.frame(lapply(densities_filt[[i]]$WT, function(x) x$x),check.names=F) %>%
      pivot_longer(cols = everything(),names_to = "Metab", values_to = "xvalue"), 
      as.data.frame(lapply(densities_filt[[i]]$WT, function(x) x$y),check.names=F) %>%
        pivot_longer(cols = everything(),names_to = "Metab", values_to = "yvalue") %>% select(yvalue))%>%
      mutate(Type = "WT"),
    # Do the same for MUT
    cbind(as.data.frame(
      lapply(densities_filt[[i]]$MUT, function(x) x$x),check.names=F) %>%
        pivot_longer(cols = everything(),names_to = "Metab", values_to = "xvalue"), 
      as.data.frame(lapply(densities_filt[[i]]$MUT, function(x) x$y),check.names=F) %>%
        pivot_longer(cols = everything(),names_to = "Metab", values_to = "yvalue") %>% select(yvalue))%>%
      mutate(Type = "MUT")) %>%
    # Add a column with the disease name
      mutate(Case = names(densities_filt)[i])
  
    d_filter = rbind(d_filter, temp)
  }
  
  cases = unique(d_filter$Case)
  case_min_max = data.frame(Case = unique(d_filter$Case))
  case_min_max$min = lapply(case_min_max$Case, function(x) min(d_filter$xvalue[d_filter$Case == x]))
  case_min_max$max = lapply(case_min_max$Case, function(x) max(d_filter$xvalue[d_filter$Case == x]))
  x_scales = list()
  # Create custom x scales for each case, regardless of the number of cases
  for(i in 1:length(cases)){
    x_scales = append(x_scales, 
           as.formula(paste0("Case == \"", cases[i], 
           "\" ~ scale_x_continuous(limits = c(", 
           unlist(case_min_max$min[case_min_max$Case == cases[i]]),",", 
           unlist(case_min_max$max[case_min_max$Case == cases[i]]),"))")))
  }
  
  metabs = d_filter %>%
    group_by(Case) %>%
    summarise("Metab" = unique(Metab))
  
  zscores_filt = semi_join(bind_rows(zscores), metabs)
  zscores_filt$colour = lapply(zscores_filt$zscore, function(x) if (x < 0) "blue" else "red")
  
  # Convert exchange reaction IDs to metabolite names for the plot
  metab_map = as.vector(metab_dict$Name)
  names(metab_map) = metab_dict$ID
  d_filter$Metab = str_replace_all(string=d_filter$Metab, pattern= fixed(metab_map))
  zscores_filt$Metab = str_replace_all(string=zscores_filt$Metab, pattern= fixed(metab_map))
  
  
  d_filter %>%
    ggplot(aes(x=xvalue,y=yvalue, fill=Type, ymin = 0, ymax = yvalue))+
    geom_ribbon(alpha=0.6)+
    facet_nested_wrap(~Case + Metab, nrow = n, scales = "free", nest_line = element_line())+
    theme_minimal()+
    ggtitle(paste0(paste0("Top flux sampling distributions filtered with a z-score threshold of ",thr)))+
    ylab("Density")+
    xlab("Flux value")+
    facetted_pos_scales(x = x_scales)+
    geom_vline(xintercept = 0, size = 0.2)+
    geom_text(data=zscores_filt, mapping = aes(col=colour),label=paste0("z-score= ",round(zscores_filt$zscore,digits = 2)), x=Inf, y=Inf, 
              inherit.aes = F, hjust = 1, vjust = 1, size = 3)
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
#' @importFrom viridis scale_colour_viridis inferno
#' @import ggrepel
#' @import ggplot2
#' @export
plot_scatter = function(sdata, zscore, means, metab_dict, thr, case="Disease"){
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
  
  metab_map = as.vector(metab_dict$Name)
  names(metab_map) = metab_dict$ID
  z_m_filtered$Metab = str_replace_all(string=z_m_filtered$Metab, pattern= fixed(metab_map))
  
  z_m_filtered = z_m_filtered %>%
    mutate(label_col = case_when(zscore>0 ~ "black",
                                 zscore<0 ~ "white"))
  
  quantile(z_m_filtered$mean_WT, probs= seq(0,1, 0.1))
  quantile(z_m_filtered$mean_MUT, probs= seq(0,1, 0.1))
  
  z_m_filtered %>% 
    ggplot(aes(x=mean_WT, y=mean_MUT, fill=zscore))+
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "white", alpha=0.5)+
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "lightgrey", alpha=0.5)+
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "lightgrey", alpha=0.5)+
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "white", alpha=0.5)+
    geom_point( size = 4, pch=21, colour = "black")+
    geom_abline(intercept = 0, linetype = "dashed")+
    # coord_cartesian(xlim = c(min(min(z_m_filtered$mean_WT),min(z_m_filtered$mean_MUT)), 
    #        max(max(z_m_filtered$mean_WT),max(z_m_filtered$mean_MUT))),
    #        ylim = c(min(min(z_m_filtered$mean_WT),min(z_m_filtered$mean_MUT)), 
    #        max(max(z_m_filtered$mean_WT),max(z_m_filtered$mean_MUT))))+
    facet_zoom(xlim = c(-1,1), ylim = c(-1,1))+
    scale_fill_gradientn(colours = c(inferno(n=7, direction = 1), inferno(n=7, direction = -1)), breaks = z_breaks_thr, labels= z_labels)+
    theme_light()+
    ggtitle(paste0("WT mean fluxes vs KO mean fluxes for ", case))+
    xlab("WT mean flux")+
    ylab("KO mean flux")+
    geom_text_repel(aes(label=Metab), colour = "black", box.padding = 0.5,min.segment.length = 0.6)
    # geom_label_repel(aes(label=Metab,fill=zscore), color= z_m_filtered$label_col)
  
}

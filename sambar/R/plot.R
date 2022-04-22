# Functions for plotting


#' Plot distributions for 1 case from raw samples
#' 
#' Compare WT and MUT flux distributions for the top set of zscore metabolites.
#' 
#' @param sdata A list containing a WT dataframe and a MUT dataframe. Each dataframe contains X samples for metabolites
#' @param zscore A dataframe of zscore calculated using calc_zscore
#' @param thr The threshold for filtering zscores
#' @param max_n The max number of metabolites to show
#' @param case The name of the disease or condition
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate 
#' @importFrom tidyr pivot_longer
#' @importFrom ggh4x facet_nested_wrap
#' @import ggplot2
#' @export
plot_distrib = function(sdata, zscore, thr, max_n=10, case="Disease"){
  #thr = 2
  #max_n = 10
  top_metab_zscore = filter_zscore(zscore, thr, max_n)
  top_metab_zscore = top_metab_zscore%>%
    mutate(colour = lapply(zscore, function(x) if (x < 0) "blue" else "red"))
  
  sdata_filter = rbind(sdata$WT %>% select(top_metab_zscore$Metab) %>% mutate(Type = "WT"),
                       sdata$MUT %>% select(top_metab_zscore$Metab) %>% mutate(Type = "MUT")) %>%
    pivot_longer(-Type, names_to = "Metab")
  
  # Change the factor levels so that the metabolites are plotted from highest zscore to lowest (absolute)
  sdata_filter$Metab <- factor(sdata_filter$Metab, levels = top_metab_zscore %>% arrange(desc(abs(zscore))) 
                               %>% pull(Metab))
  
  gp = sdata_filter %>% mutate(Case = case) %>%
    ggplot(aes(x=value, fill=Type))+
    geom_density(alpha = 0.5, colour=NA)+
    theme_minimal()+
    geom_vline(xintercept = 0, size = 0.2)+
    facet_nested_wrap(Case ~ factor(Metab), nrow =1, scales = "free_y", strip.position = "top", 
                      nest_line = element_line())+
    geom_text(data=top_metab_zscore, mapping = aes(col=colour),
              label=paste0("z-score= ",round(top_metab_zscore$zscore,digits = 2)), x=Inf, y=Inf,
              inherit.aes = F, hjust = 1, vjust = 1, size = 3)+
    ggtitle(paste0(paste0("Top ",max_n," flux sampling distributions filtered with a z-score threshold of ",thr,
                          " for ",case)))+
    ylab("Density")+
    xlab("Flux value")
  
  return(gp)
}

#' Plot a basic demo distribution
#' 
#' Compare WT and MUT flux distributions for the top set of zscore metabolites for an example disease.
#' 
#' @export
plot_demo_distrib = function(){
  # Generate hypothetical data
  example_d = list(
    "WT"=data.frame("MET_X"=rnorm(5000, 0, 1)), 
    "MUT"=data.frame("MET_X"=rnorm(5000, 4, 1)))
  zscore = data.frame(Metab="MET_X", zscore=2.5)
  
  plot_distrib(example_d, zscore,thr = 2,case = "Example")
}



#' Plot multiple distributions from raw samples
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


#' Plot multiple density distributions for different diseases
#' 
#' Compare WT and MUT flux distributions for the top set of zscore metabolites for multiple diseases using filtered 
#' densities as inputs instead of raw sdata.
#' 
#' @param densities_filt A list of case lists, each containing a WT list and a MUT list. Each WT and MUT list contains
#'  X densities for metabolites
#' @param metab_dict A dataframe containing an ID column and a name column. This is useful for converting metabolite IDs
#' to readable names. It can be generated using sambapy from the metabolic model.
#' @param thr A numeric threshold. This will be used in the title for plotting purposes only.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select
#' @importFrom tidyr pivot_longer
#' @import grid
#' @import ggplot2
#' @import gtable
#' @importFrom ggh4x facet_nested_wrap facetted_pos_scales
#' @export
plot_multi_density_distrib = function(densities_filt, metab_dict, thr=2){
  #densities_filt = d_all_filtered
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
  # Join the metabs with zscores and add colour for the plot
  metabs = d_filter %>%
    group_by(Case) %>%
    summarise("Metab" = unique(Metab))
  
  zscores_filt=zscores %>% bind_rows(.id = "Case") %>%
    semi_join(.,metabs, by=c("Case", "Metab")) %>%
    mutate(colour = lapply(zscore, function(x) if (x < 0) "blue" else "red"))
  
  # Convert exchange reaction IDs to metabolite names for the plot
  metab_map = as.vector(metab_dict$Name)
  names(metab_map) = metab_dict$ID
  d_filter$Metab = str_replace_all(string=d_filter$Metab, pattern= fixed(metab_map))
  zscores_filt$Metab = str_replace_all(string=zscores_filt$Metab, pattern= fixed(metab_map))
  
  # Get number of cases for n. of rows for the plot
  n_cases = length(cases)
  n_rows = n_cases + 2
  # Create gtable object
  g = gtable(widths = unit(c(1,6,1), c("lines", "null", "null")),
             heights = unit(c(1,rep(1,n_cases),1), c("line",rep("null",n_cases), "line")))
  # Add titles
  title = paste0("Top flux sampling distributions filtered with a z-score threshold of ",thr)
  g = gtable_add_grob(g, textGrob(title), t=1,b=1,l=2,r=2)
  g = gtable_add_grob(g, textGrob("Flux value"), t=n_rows ,b=n_rows ,l=2,r=2)
  g = gtable_add_grob(g, textGrob("Density", rot=90), t=1,b=n_rows ,l=1,r=1)
  
  # Create a temp plot to extract the legend
  gp = d_filter %>% filter(Case == cases[1]) %>%
    ggplot(aes(x=xvalue,y=yvalue, fill=Type, ymin = 0, ymax = yvalue))+
    geom_ribbon(alpha=0.6)
  # Extract the legend
  guide = gtable_filter(ggplotGrob(gp), pattern="guide")
  # Add the legend to the table object
  g = gtable_add_grob(g,guide , t=1,b=n_rows,l=3,r=3)
  
  # Loop over all cases
  for(c in 1:length(cases)){
    # Filter for the current case
    zscores_c = zscores_filt %>% filter(Case == cases[c])
    d_filter_c = d_filter %>% filter(Case == cases[c])
    # Change the factor levels so that the metabolites are plotted from highest zscore to lowest (absolute)
    d_filter_c$Metab <- factor(d_filter_c$Metab, levels = zscores_c %>% arrange(desc(abs(zscore))) %>% pull(Metab))
    
    # Create a plot for the current case
    gp = d_filter_c %>%
      ggplot(aes(x=xvalue,y=yvalue, fill=Type, ymin = 0, ymax = yvalue))+
      geom_ribbon(alpha=0.6)+
      facet_nested_wrap(Case ~ factor(Metab), nrow =1, scales = "free_y", strip.position = "top", nest_line = element_line())+
      theme_minimal()+
      geom_vline(xintercept = 0, size = 0.2)+
      geom_text(data=zscores_c, mapping = aes(col=colour),
                label=paste0("z-score= ",round(zscores_c$zscore,digits = 2)), x=Inf, y=Inf,
                inherit.aes = F, hjust = 1, vjust = 1, size = 3)+ 
      # Remove titles, axis labels and legend since we've already added them in the main plot
      labs(title = element_blank())+ 
      theme(axis.title.x = element_blank())+ 
      theme(axis.title.y = element_blank())+ 
      theme(legend.position="none")
    
    # Add the plot panel to table object
    g <- gtable_add_grob(g,ggplotGrob(gp), t=(c+1),b=(c+1),l=2,r=2)
    
  }
  
  # Plot
  grid.newpage()
  grid.draw(g)
}


#' Plot a scatter plot of WT means vs MUT means
#' 
#' Plots a scatter plot with WT flux means on the x-axis and MUT flux means on the y-axis. The points are coloured
#' by zscore.
#' 
#' @param sdata A list containing a WT dataframe and a MUT dataframe. Each dataframe contains X samples for metabolites
#' @param zscore A dataframe of zscore calculated using calc_zscore
#' @param means A dataframe of means calculated using calc_sdata_means
#' @param metab_dict A dataframe containing an ID column and a name column. This is useful for converting metabolite IDs
#' to readable names. It can be generated using sambapy from the metabolic model.
#' @param thr The threshold for filtering zscores
#' @param case The name of the disease or condition
#' @param outlier.dist The distance over which a metabolite is considered an outlier for plotting purposes.
#' @param outlier.count The minimum number of outliers there needs to be for a plot zoom to be added.
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom tidyr gather pivot_longer
#' @importFrom ggforce facet_zoom
#' @importFrom tibble rownames_to_column
#' @importFrom viridis scale_colour_viridis inferno
#' @import ggrepel
#' @import ggplot2
#' @export
plot_scatter = function(sdata, zscore, means, metab_dict, thr, case="Disease", outlier.dist=5, outlier.count=0){
  # Join the zscore with the means
  z_m = zscore %>% left_join(., means$WT , by="Metab") %>% left_join(., means$MUT , by="Metab")
  # Filter the zscore based on the provided threshold
  z_m_filtered = z_m %>% filter(zscore > thr | zscore < -thr)
  # Create the breaks for the legend
  z_max = max(abs(z_m_filtered$zscore))
  z_breaks = round(seq(from=-z_max, to=z_max,by=z_max/2))
  z_breaks = z_breaks[abs(z_breaks-thr) > 0.9]
  z_breaks_thr = c(z_breaks, thr, -thr)
  z_labels = c(z_breaks, paste0(" thr=",thr), paste0("-thr=", -thr))
  # Create the scale limits for the legend so that the scale gradient is symmetrical
  z_limits = c(min(z_breaks),max(z_breaks))
  
  metab_map = as.vector(metab_dict$Name)
  names(metab_map) = metab_dict$ID
  z_m_filtered$Metab = str_replace_all(string=z_m_filtered$Metab, pattern= fixed(metab_map))
  
  z_m_filtered = z_m_filtered %>%
    mutate(label_col = case_when(zscore>0 ~ "black",
                                 zscore<0 ~ "white"))
  # Check for outliers
  all_zscores = c(z_m_filtered$mean_WT, z_m_filtered$mean_MUT)
  q_z = quantile(all_zscores)
  iqr_z = IQR(all_zscores)
  # More relaxed bounds: used for plotting
  z_thr_upper_med = (iqr_z * 1.5) + q_z[4]
  z_thr_lower_med = q_z[2] - (iqr_z * 1.5)
  # Stricter bounds: used for the actual check
  z_thr_upper = (iqr_z * 3) + q_z[4]
  z_thr_lower = q_z[2] - (iqr_z * 3)
  # If there are any values further than 5 from the upper quantile check, add a zoom
  if (sum(unlist(lapply(all_zscores, function(x) abs(x-z_thr_upper))) > outlier.dist) > outlier.count){
    zoom = T
  }else{
    zoom = F
  }
  
  
  z_m_filtered %>% 
    ggplot(aes(x=mean_WT, y=mean_MUT, fill=zscore))+
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "white", alpha=0.5)+
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "lightgrey", alpha=0.5)+
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "lightgrey", alpha=0.5)+
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "white", alpha=0.5)+
    geom_point( size = 4, pch=21, colour = "black")+
    geom_abline(intercept = 0, linetype = "dashed")+ {
      if (zoom == T){
        facet_zoom(xlim = c(z_thr_lower_med,z_thr_upper_med), ylim = c(z_thr_lower_med,z_thr_upper_med))
      }
    }+
    scale_fill_gradientn(colours = c(inferno(n=7, direction = 1), inferno(n=7, direction = -1)), 
                         breaks = z_breaks_thr, labels= z_labels,
                         limits = z_limits)+
    theme_light()+
    ggtitle(paste0("WT mean fluxes vs KO mean fluxes for ", case))+
    xlab("WT mean flux")+
    ylab("KO mean flux")+
    geom_text_repel(aes(label=Metab), colour = "black", box.padding = 0.5,min.segment.length = 0.6)
  
}

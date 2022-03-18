# library(stringr)
library(tidyverse)
library(sambar)
library(ggridges)
library(ggforce)
library(viridis)
library(tools)
#pred = data.frame(Metab=metab_names$Metab.EX.IDs)
name = "wtopt_0_1"

ko = data.table::fread(paste0("/home/juliette/these/data/IEM/sampling/results/",name,"_100000/ko/xaa_Recon-2_from_matlab_sampling_100000_KO.csv.gz"))
wt = data.table::fread(paste0("/home/juliette/these/data/IEM/sampling/results/",name,"_100000/wt/xaa_Recon-2_from_matlab_sampling_100000_base.csv.gz"))
sdata = list(ko,wt)
names(sdata) = c("c2","c1")
sdata = lapply(sdata, as.data.frame)

diff = calc_diff(sdata)
#write.csv(diff, "/home/juliette/these/data/samba_tests/results/diff/xaa_diff.csv")
#diff =  data.table::fread("/home/juliette/these/data/samba_tests/results/diff/xaa_diff.csv",check.names = F)
# diff = diff[,2:3]
#write.csv(zscore, "/home/juliette/these/data/samba_tests/results/zscore/xaa_zscore.csv")
zscore = read.csv("/home/juliette/these/data/samba_tests/results/zscore/xaa_zscore.csv",check.names = F)
zscore = zscore[,2:3]

# Scatter plot
avg_c1 = sdata$c1 %>%
  summarise(across(everything(), mean))
avg_c2 = sdata$c2 %>%
  summarise(across(everything(), mean))
t1 = avg_c1 %>% t() %>% as.data.frame() %>%
  rownames_to_column("Metab") %>%
  rename("mean_c1" = V1)
t2 = avg_c2 %>% t() %>% as.data.frame() %>%
  rownames_to_column("Metab") %>%
  rename("mean_c2" = V1)
z_m = zscore %>% left_join(., t1 , by="Metab") %>% left_join(., t2 , by="Metab")

thr = 1
z_m_filtered = z_m %>% filter(zscore > thr | zscore < -thr)
z_breaks = round(seq(from=min(z_m_filtered$zscore), to=max(z_m_filtered$zscore),by=5))
z_breaks = z_breaks[abs(z_breaks-thr) > 4]
z_breaks_thr = c(z_breaks, thr, -thr)
z_labels = c(z_breaks, paste0(" thr=",thr), paste0("-thr=", -thr))

z_m_filtered %>% 
  ggplot(aes(x=mean_c1, y=mean_c2, colour=zscore))+
  geom_point(size = 4, alpha = 0.8)+
  geom_abline(intercept = 0, linetype = "dashed")+
  # coord_cartesian(xlim = c(min(min(z_m_filtered$mean_c1),min(z_m_filtered$mean_c2)), 
  #        max(max(z_m_filtered$mean_c1),max(z_m_filtered$mean_c2))),
  #        ylim = c(min(min(z_m_filtered$mean_c1),min(z_m_filtered$mean_c2)), 
  #        max(max(z_m_filtered$mean_c1),max(z_m_filtered$mean_c2))))+
  facet_zoom(xlim = c(-1,1), ylim = c(-1,1))+
  scale_colour_viridis(option = "magma", breaks = z_breaks_thr, labels= z_labels)+
  theme_light()+
  ggtitle("WT mean fluxes vs KO mean fluxes")+
  xlab("WT mean flux")+
  ylab("KO mean flux")




sd(zscore$zscore,na.rm = T)
hist(zscore$zscore, breaks = 50)
thr = 0.01

non_signif = zscore[zscore$zscore < thr & zscore$zscore > -thr,]



thr = 2
max_n = 10
top_metab_zscore = zscore %>% arrange(desc(abs(zscore))) %>%
  slice_head(n = max_n) %>% filter(abs(zscore) > thr)

ko_filter = ko %>% select(top_metab_zscore$Metab)
wt_filter = wt %>% select(top_metab_zscore$Metab)

sdata_filter = rbind(ko_filter %>% mutate(Type = "c2"),
                    wt_filter %>% mutate(Type = "c1")) %>%
  pivot_longer(-Type, names_to = "Metab")

sdata_filter %>% 
  ggplot(aes(x=value, y=factor(Metab, levels = rev(top_metab_zscore$Metab)), fill=Type, height= ..ndensity..))+
  geom_density_ridges(rel_min_height = 0.01, alpha=0.5, scale=1.5, linetype="blank")+
  scale_fill_manual(values = c("firebrick1", "dodgerblue"))+
  coord_cartesian(clip = "off")+
  theme_ridges(center_axis_labels = T)+
  ggtitle(paste0("Flux sampling distributions"))+
  ylab("Metabolites")+
  xlab("Flux value")

# Generic example distrib
example_d = data.frame(c1=rnorm(5000, 0, 1), c2=rnorm(5000, 2, 1))
example_d %>% pivot_longer(cols = c(c1,c2),names_to = "Type", values_to = "value") %>%
  ggplot(aes(x=value, fill=Type))+
  geom_density(alpha = 0.5, colour=NA)+
  theme_minimal()+
  ggtitle(paste0("Flux sampling distributions"))+
  ylab("Disease")+
  xlab("Flux value")
  

# Flux distrib on one line
sdata_filter %>%
  ggplot(aes(x=value, fill=Type))+
  geom_density(alpha = 0.5, colour=NA)+
  facet_wrap(~Metab,nrow = 1, scales = "free")+
  theme_minimal()+
  ggtitle(paste0("Flux sampling distributions"))+
  ylab("Disease")+
  xlab("Flux value")


zscore_art = readxl::read_excel("~/Downloads/journal.pcbi.1009522.s006(1).xlsx",skip = 1)
hist(zscore_art$`Z-score`, breaks = 50)


# Compare variability of zscore WT vs WT and WT vs KO
wt_1 = sdata$c1[sample(nrow(sdata$c1), 25000),]
wt_2 = sdata$c1[sample(nrow(sdata$c1), 25000),]
sdata = list(wt_2, wt_1)
names(sdata) = c("c2","c1")
sdata = lapply(sdata, as.data.frame)
(sum(sdata$c1[,1] <= sdata$c2[,1])+1)/(25000+1)
(sum(sdata$c1[,1] > sdata$c2[,1])+1)/(25000+1)

sd(sdata$c1[,1])
sd(sdata$c2[,1])

diff_wt = calc_diff(sdata)
zscore_wt = calc_zscore(diff_wt)
hist(zscore_wt$zscore, breaks = 50)
sd(zscore_wt$zscore, na.rm = T)
sd(zscore$zscore, na.rm = T)

library(coin)
independence_test(sdata$c1[,1] ~ sdata$c2[,1])
independence_test(sdata$c1[,1] ~ sdata$c2[,1])


shapiro.test(zscore_wt$zscore)
# full gaussienne
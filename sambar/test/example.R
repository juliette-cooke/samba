library(sambar)
library(ggh4x)
library(cowplot)
library(bigmemory)
library(tidyverse)
# To import one disease (two files)
sdata = load_sampling_results("/home/juliette/these/code/git/samba/data/", "xaa")

diff = calc_diff(sdata)
zscore = calc_zscore(diff)
gc()

plot_distrib(sdata, zscore, thr=2, max_n = 10)
plot_scatter(sdata, zscore, means, thr=1)


# To plot more than one disease in 1 figure
# Put all disease files in the same folder
# And add prefixes (=disease names at the beginning of your files)
path_to_folder = "/home/juliette/these/code/git/samba/data/"
prefixes = c("xaa", "xab")
zscores = list()
sdata_all_filtered = list()

for (i in 1:length(prefixes)){
  sdata = load_sampling_results(path_to_folder, prefixes[i])
  diff = calc_diff(sdata, divfactor = 1)
  gc()
  zscores[[i]] = calc_zscore(diff)
  names(zscores)[i] = prefixes[i]

  # Select filtered metabs
  thr = 2
  max_n = 10
  metabs_filtered = zscores[[i]] %>% arrange(desc(abs(zscore))) %>%
    slice_head(n = max_n) %>% filter(abs(zscore) > thr)
  sdata_all_filtered[[i]] = lapply(sdata, function(x) x %>% select(metabs_filtered$Metab))
  names(sdata_all_filtered)[i] = prefixes[i]
}

demo = plot_basic_distrib()
distribs = plot_multi_distrib(sdata_all_filtered, 2)

plot_grid(plotlist = list(demo, distribs), nrow = 2,scale=c(0.5,1), align = "v",axis = "t",rel_heights = c(0.25,1))



# Testing density
sdata = load_sampling_results("/home/juliette/these/code/git/samba/data/", "xaa")
# hist(sdata$MUT$`EX_10fthf(e)`, breaks = 100)
t = density(sdata$MUT$`EX_10fthf(e)`)

# Plotting samples without plotting all 100 000 samples (slow)
ggplot(NULL, aes(x=t$x, y=t$y))+
  geom_line()


# random.points = approx(
#   cumsum(t$y)/sum(t$y),
#   t$x,
#   runif(10000))$y
# hist(random.points, breaks = 100)
path_to_folder = "/home/juliette/these/code/git/samba/data/"
prefixes = c("xaa", "xab")
zscores = list()
d_all_filtered = list()

i=1
for (i in 1:length(prefixes)){
  sdata = load_sampling_results(path_to_folder, prefixes[i])
  diff = calc_diff(sdata)
  gc()
  zscores[[i]] = calc_zscore(diff)
  names(zscores)[i] = prefixes[i]
  
  # Select filtered metabs
  thr = 2
  max_n = 10
  metabs_filtered = zscores[[i]] %>% arrange(desc(abs(zscore))) %>%
    slice_head(n = max_n) %>% filter(abs(zscore) > thr)
  densities = calc_density(sdata)
  d_all_filtered[[i]] = lapply(densities, function(x) x)
  d_all_filtered[[i]] = list(densities$WT[names(densities$WT) %in% metabs_filtered$Metab],
                        densities$MUT[names(densities$MUT) %in% metabs_filtered$Metab])
  names(d_all_filtered)[i] = prefixes[i]
  names(d_all_filtered[i][[prefixes[i]]]) = c("WT", "MUT")
}

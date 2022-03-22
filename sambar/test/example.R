library(sambar)
library(cowplot)
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
  diff = calc_diff(sdata)
  gc()
  zscores[[i]] = calc_zscore(diff)
  names(zscores)[i] = prefixes[i]

  # Select filtered metabs
  sdata_all_filtered[[i]] = calc_filter_sdata(sdata, zscores[[i]], thr = 2, max_n = 10)
  names(sdata_all_filtered)[i] = prefixes[i]
}

# Add demo plot
demo = plot_basic_distrib()
distribs = plot_multi_distrib(sdata_all_filtered, 2)
plot_grid(plotlist = list(demo, distribs), nrow = 2,scale=c(0.5,1), align = "v",axis = "t",rel_heights = c(0.25,1))


# Plot multiple density plots
path_to_folder = "/home/juliette/these/code/git/samba/data/"
prefixes = c("xaa", "xab")
zscores = list()
d_all_filtered = list()

for (i in 1:length(prefixes)){
  sdata = load_sampling_results(path_to_folder, prefixes[i])
  diff = calc_diff(sdata)
  gc()
  zscores[[i]] = calc_zscore(diff)
  names(zscores)[i] = prefixes[i]
  
  # Select filtered metabs
  d_all_filtered[[i]] = calc_filter_density(sdata, zscores[[i]], thr = 2, max_n = 10)
  names(d_all_filtered)[i] = prefixes[i]
  names(d_all_filtered[i][[prefixes[i]]]) = c("WT", "MUT")
}

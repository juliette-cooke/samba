library(sambar)
library(cowplot)
library(tidyverse)
library(tools)

# Read in metab dict for the specific model
metab_dict = read.csv("/home/juliette/these/code/git/samba/test_data/Recon-2_from_matlab_metab_id_name.tsv", sep = "\t")


# To import one disease (two files)
sdata = load_sampling_results("/home/juliette/these/code/git/samba/data/", "xac")

diff = calc_diff(sdata)
zscore = calc_zscore(diff)
gc()
plot_distrib(sdata, zscore, thr=2, max_n = 10)
means = calc_sdata_means(sdata)
plot_scatter(sdata, zscore, metab_dict = metab_dict, means, thr=1,outlier.dist = 3)


# Create a loop producing scatter plots for multiple diseases
# Automatically parse the prefixes from files in the provided folder
indir="/home/juliette/these/code/git/samba/data/"
files = list.files(indir)
prefixes = unlist(unique(lapply(files, function(x) unlist(str_split(file_path_sans_ext(x, compression = T), "_"))[1])))
prefixes = prefixes[1:3]
for (i in 1:length(prefixes)){
  sdata = load_sampling_results("/home/juliette/these/code/git/samba/data/", prefixes[i])
  diff = calc_diff(sdata)
  zscore = calc_zscore(diff)
  gc()
  means = calc_sdata_means(sdata)
  plot_scatter(sdata, zscore, metab_dict = metab_dict, means, thr=1)
  ggsave(paste0("/home/juliette/these/code/git/samba/test_plots/scatter_all/",prefixes[i], "_scatterplot.png"),width = 15, height = 10)
}


# To plot more than one disease in 1 figure
# This can take a while with multiple diseases as we're storing 100 000 samples x 2 for each disease
# If you have many diseases it's better to calculate density estimations before plotting (see below)
# Put all disease files in the same folder
# And add prefixes (=disease names/IDs at the beginning of your files)
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



# Plot multiple density plots: more efficient to calculate density approximations for each sampling
# Read in metab dict for the specific model
metab_dict = read.csv("/home/juliette/these/code/git/samba/test_data/Recon-2_from_matlab_metab_id_name.tsv", sep = "\t")
# Read in and convert each disease into density plots (uses less RAM and plots easier)
path_to_folder = "/home/juliette/these/code/git/samba/data/"
files = list.files(path_to_folder)
prefixes = unlist(unique(lapply(files, function(x) unlist(str_split(file_path_sans_ext(x, compression = T), "_"))[1])))
prefixes = prefixes[1:3]
zscores = list()
d_all = list()
for (i in 1:length(prefixes)){
  sdata = load_sampling_results(path_to_folder, prefixes[i])
  diff = calc_diff(sdata)
  gc()
  zscores[[i]] = calc_zscore(diff)
  names(zscores)[i] = prefixes[i]
  
  # Calculate densities for each metabolite (no filtering)
  d_all[[i]] = calc_density(sdata)
  names(d_all)[i] = prefixes[i]
}

# Filter outside the loop
thr = 2
d_all_filtered = filter_all_density(d_all, zscores, thr = thr, max_n = 10)

plot_multi_density_distrib(d_all_filtered, metab_dict, thr = thr)

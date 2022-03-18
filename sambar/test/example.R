library(sambar)

# To import one disease (two files)
sdata = load_sampling_results("/home/juliette/these/code/git/samba/data/", "xaa")

diff = calc_diff(sdata, divfactor = 1)
zscore = calc_zscore(diff)
gc()
plot_distrib(sdata, zscore, thr=2, max_n = 10)

# To import more than one disease into 1 figure


plot_scatter(sdata, zscore, thr=1)

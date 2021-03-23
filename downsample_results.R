library(glue)
library(tidyr)
library(ggplot2)
library(dplyr)

my_cell_types = c("OPC","Astrocytes","EN-MIX","EN-PYR","Glia-MIX","IN-SST","IN-SV2C",
                  "IN-VIP","Microglia","OL","IN-PVALB",
                  "EN-L2-3", "EN-L4", "EN-L5-6")

my_sims = 1:25
out_dir = "/home/jgagnon1/scRNAseq/downsample"
prefix = "MS_Nature_downsample_0R1"

DEG_table = sapply(my_cell_types,function(cluster)
{
  print(cluster)
  sapply(my_sims,function(x) {
    my_data = read.csv(file.path(out_dir,glue("{prefix}_glmmTMB_cell_{cluster}_wo_cdr_sim_{x}.csv")))
    c(sum(my_data$FDR < 0.05, na.rm = TRUE))
   } 
  )
  }
)

plot_data = data.frame(DEG_table) %>% gather(cell_type,nDEG)
mean_data = plot_data %>% group_by(cell_type) %>% summarise(my_mean=mean(nDEG),my_sd=sd(nDEG))

ggplot(mean_data, aes(x=reorder(cell_type,-my_mean), y=my_mean, fill=cell_type)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=pmax(my_mean-2*my_sd,0), ymax=my_mean+2*my_sd), width=.2)+
  xlab("cell type")+ylab("Average DEGs +/- 2 sd")+ggtitle("glmmTMB: MS vs Control with downsampling")

pDEG_table = sapply(my_cell_types,function(cluster)
{
  print(cluster)
  sapply(my_sims,function(x) {
    my_data = read.csv(file.path(out_dir,glue("{prefix}_glmmTMB_cell_{cluster}_wo_cdr_sim_{x}.csv")))
    c(sum(my_data$FDR < 0.05, na.rm = TRUE)/nrow(my_data)*100)
  } 
  )
}
)

plot_data = data.frame(pDEG_table) %>% gather(cell_type,pDEG)
mean_data = plot_data %>% group_by(cell_type) %>% summarise(my_mean=mean(pDEG),my_sd=sd(pDEG))

ggplot(mean_data, aes(x=reorder(cell_type,-my_mean), y=my_mean, fill=cell_type)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=pmax(my_mean-2*my_sd,0), ymax=my_mean+2*my_sd), width=.2)+
  xlab("cell type")+ylab("Average %DEGs +/- 2 sd")+ggtitle("glmmTMB: MS vs Control with downsampling")

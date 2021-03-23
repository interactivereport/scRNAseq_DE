library(ggplot2)
library(dplyr)

cell_type = "Astrocytes"
data_dir= "/home/jgagnon1/scRNAseq/300k_results_CALvsNWM/"

my_results = lapply(seq(16), function(x) read.csv(paste0(data_dir,"NWMCAL_glmmTMB_cell_",cell_type,"_",x,".csv")))
merged = do.call(rbind,my_results)
merged$FDR = p.adjust(merged$Pvalue,"fdr")
merged$sig = ifelse(merged$FDR <0.05,"sig","ns")

ngenes = nrow(merged)
nDE = sum(merged$FDR < 0.05,na.rm=T)
ngenes
nDE
nDE/ngenes*100 

merged %>% ggplot(aes(log2FC,-log10(Pvalue)))+geom_point(aes(color=sig))+theme_bw()+ggtitle(cell_type)
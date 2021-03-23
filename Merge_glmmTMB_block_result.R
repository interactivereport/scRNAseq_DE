
library(ggplot2)
library(dplyr)
library(knitr)
library(kableExtra)
# block results path
data_dir = "/home/qwan/scRNAseq/RESULTS/DEV/350kDE_REPORT/NAWMsubsetSPMSvsPPMS_meta_filtered/"  

# set up cell types
my_cell_types = c("Astrocytes","EN-L4","EN-L2-3", "EN-L5-6","EN-MIX","EN-PYR","Glia-MIX","IN-SST","IN-SV2C","IN-VIP","Microglia","OL","OPC","IN-PV")
ntypes = length(my_cell_types)

# create a data set to save the summary table
my_frame = data.frame("nDE"=rep(NA,ntypes),"%DE"=rep(NA,ntypes),"%fail"=rep(NA,ntypes),"ngenes"=rep(NA,ntypes),check.names = FALSE,row.names = my_cell_types)

for (cell_type in my_cell_types){
  my_results = lapply(seq(16), function(x) read.csv(paste0(data_dir,"SPMSPPMS_glmmTMB_cell_",cell_type,"_",x,".csv"))) 
  merged = do.call(rbind,my_results)
  merged$FDR = p.adjust(merged$Pvalue,"fdr") #important
  merged$sig = ifelse(merged$FDR <0.05,"sig","ns")
  
  ngenes = nrow(merged)
  nDE = sum(merged$FDR < 0.05,na.rm=T)
  pDE = round(nDE/ngenes*100,1)
  pFail = round(sum(merged$Error)/nrow(merged)*100,1)
  
  my_frame[cell_type,"nDE"] = nDE
  my_frame[cell_type,"%DE"] = pDE
  my_frame[cell_type,"%fail"] = pFail
  my_frame[cell_type,"ngenes"] = ngenes
  
  # save merge result as csv file
  write.csv(merged,
            file = file.path("/home/qwan/scRNAseq/RESULTS/DEV/350kDE_REPORT/CALvsNWMdiedAge/", glue("New_350k_glmmTMB_cell_{cell_type}.csv")),
            row.names = F)
  
  #print(merged %>% dplyr:::filter(!is.na(sig)) %>% ggplot(aes(log2FC,-log10(Pvalue)))+geom_point(aes(color=sig))+theme_bw()+labs(title=cell_type,subtitle=paste0("%DE: ", pDE,"% DE: ", nDE," Genes: ", ngenes, " Failed: ", pFail,"%")))
  #print(merged %>% dplyr:::filter(FDR < 0.05) %>% dplyr:::select(ID,log2FC,Pvalue,FDR) %>% dplyr:::arrange(desc(abs(log2FC))) %>% dplyr:::mutate(FC=ifelse(log2FC > 0, 2^log2FC,2^(-log2FC)))%>% head(20))
  
}

my_frame<-cbind(my_frame,my_frame1)

kable(my_frame) %>%  kable_styling(bootstrap_options = "striped", full_width = F)
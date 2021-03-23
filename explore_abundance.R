suppressMessages(library(tidyverse))

source("pipeline_class_042420.R")

opt = list(in_dir = "/home/jgagnon1/jupyter-notebook-dir", gene_info_file = "MS_Nature_gene_info.csv", 
           meta_file = "MS_Nature_meta_data.csv", sampleID_var = "sample", cluster_var = "cell_type",
           group_var = "diagnosis", reference_group = "Control", alternative_group = "MS", 
           count_file = "MS_Nature_UMI_subset.mtx", covars = "Capbatch,age,sex", 
           out_dir =  "/home/jgagnon1/scRNAseq/", 
           out_prefix = "MS_Nature",glmm_model="nbinom2")

if (is.na(opt$covars)) {
  covars = NULL 
} else {
  covars = unlist(strsplit(opt$covars, split = "[ ,;]+"))
}

UMI_data<-readMM(file.path(opt$in_dir, opt$count_file))
meta_info=read.csv(file.path(opt$in_dir, opt$meta_file), stringsAsFactors=FALSE)
gene_info=read.csv(file.path(opt$in_dir, opt$gene_info_file), stringsAsFactors=FALSE)
meta_info<-meta_info%>%dplyr:::mutate(cell_type= case_when(cell_type=="OL-C"|cell_type=="OL-B"|cell_type=="OL-A"~"OL",
                                                           cell_type=="EN-L2-3-A"|cell_type=="EN-L2-3-B"~"EN-L2-3",
                                                           TRUE~cell_type))
meta_info$cell_type <- gsub(" ", "_", meta_info$cell_type)
meta_info<-meta_info%>%dplyr:::mutate(stage= case_when(stage=="Acute/Chronic active"~"Active",
                                                       stage=="Chronic inactive"~"Inactive",
                                                       TRUE~"Control"))
count_data <- t(UMI_data)
colnames(count_data)<-meta_info$cell
rownames(count_data) <- gene_info$index

# Load the data
# load the data
sce <- BiostatsSingleCell$new(count_data = count_data,
                              meta_data = meta_info,
                              sampleId_col = opt$sampleID_var,
                              cluster_col = opt$cluster_var,
                              treatment_col = opt$group_var)

# 1st round of filtering
sce$apply_filter(min.perc.cells.per.gene = 0.00) # 0% expression requirement

all_cell_types = unique(sce$pData()[,"cell_type"])
all_cell_types = na.omit(all_cell_types)
print(all_cell_types)

for (cell_type in all_cell_types)
{
  print(cell_type)
  error_flag <- FALSE
  
  tryCatch({
    sce$set_group_mode(cell_type,opt$reference_group,opt$alternative_group)
    }, error = function(err) {
    error_flag <<- TRUE
    }
  )
  if (error_flag) next
  
  tryCatch({
    filterR2 = sce$apply_filter_contrasts_R6()
    }, error = function(err) {
    error_flag <<- TRUE
    }
  )
  if (error_flag) next
  print(dim(filterR2$assayData()))
}
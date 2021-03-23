suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

source("pipeline_class_042420.R")
option_list = list(
  make_option("--in_dir", action = "store", default = NA, type = "character",
              help = "Path to input files [required]"),
  make_option("--gene_info_file", action = "store", default = NA, type = "character",
              help = "File name of the gene info file [required]"),
  make_option("--meta_file", action = "store", default = NA, type = "character",
              help = "File name of the meta file [required]"),
  make_option("--sampleID_var", action = "store", default = "sample", type = "character",
              help = "The name of the sample ID variable in the meta file [default: %default]"),
  make_option("--cluster_var", action = "store", default = "cell_type", type = "character",
              help = "The name of the cluster variable in the meta file [default: %default]"),
  make_option("--cluster", action = "store", default = "cell_type", type = "character",
              help = "Which cluster to use? [default: %default]"),
  make_option("--group_var", action = "store", default = "diagnosis", type = "character",
              help = "The name of the group variable in the meta file [default: %default]"),
  make_option("--reference_group", action = "store", default = "Control", type = "character",
              help = "The name of the reference group [default: %default]"),
  make_option("--minimum_percentage_cell", action = "store", default = 0.1, type = "double",
              help = "minimum percentage cell expressed [default: %default]"),
  make_option("--alternative_group", action = "store", default = "MS", type = "character",
              help = "The name of the alternative group [default: %default]"),
  make_option("--covars", action = "store", default = NA, type = "character",
              help = "The names of covars to be included for the analysis"),
  make_option("--covars_formula", action = "store", default = NA, type = "character",
              help = "The formula of covars to be included for the analysis"),
  make_option("--count_file", action = "store", default = NA, type = "character",
              help = "File name of the count data [required]"),
  make_option("--out_dir", action = "store", default = NA, type = "character",
              help = "Path for output file [required]"),
  make_option("--out_prefix", action = "store", default = NA, type = "character",
              help = "Prefix of the output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)

# opt = list(in_dir = "/home/jgagnon1/jupyter-notebook-dir", gene_info_file = "MS_Nature_gene_info.csv", 
#            meta_file = "MS_Nature_meta_data.csv", sampleID_var = "sample", cluster_var = "cell_type",
#            cluster = "Astrocytes", group_var = "diagnosis", reference_group = "Control", alternative_group = "MS", 
#            count_file = "MS_Nature_UMI_subset.mtx", covars = "Capbatch,age,sex", 
#            out_dir =  "/home/jgagnon1/scRNAseq/downsample/", 
#            out_prefix = "MS_Nature",glmm_model="nbinom2")

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
sce$apply_filter(min.perc.cells.per.gene = 0.00) 

# set mode
sce$set_group_mode(cluster_of_interest = opt$cluster, ref_group = opt$reference_group, alt_group = opt$alternative_group)

# Filtering round 2
sce_qc <- sce$apply_filter_contrasts_R6()

cluster = opt$cluster
prefix = opt$out_prefix

# 25 downsamples
nsim = 25

for(sim in seq(nsim))
{
  start_time = proc.time()
  sce_down <- sce_qc$down_sample(775)
  glmmTMB.res <- sce_down$glmmTMB_pipeline(covs = covars, family = "nbinom2", detection_rate = FALSE,cores=8)
  cat("nDEGs: ",sum(glmmTMB.res$FDR < 0.05,na.rm=TRUE),"\n")
  write.csv(glmmTMB.res, file.path(opt$out_dir, glue("{prefix}_glmmTMB_cell_{cluster}_wo_cdr_sim_{sim}.csv")), row.names = FALSE)
  print(proc.time() - start_time)
}

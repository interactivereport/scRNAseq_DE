# ------------------------------------------------------------------------------------
# Set up library path & load libraries
# ------------------------------------------------------------------------------------

#.libPaths(c("/camhpc/home/lhou1/R/x86_64-pc-linux-gnu-library/3.5","/camhpc/pkg/R/3.5.0/centos7/lib64/R/library"))
# .libPaths(c("/camhpc/home/lhou1/R/x86_64-pc-linux-gnu-library/3.5",
#             "/camhpc/pkg/R/3.5.1/centos6/lib64/R/library",
#             "/camhpc/home/qwan/R/x86_64-pc-linux-gnu-library/3.5"
#             ))

suppressMessages(library(data.table))
suppressMessages(library(glue))
suppressMessages(library(scater))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(Matrix))
suppressMessages(library(MAST))
suppressMessages(library(fst))
suppressMessages(library(tidyr))
suppressMessages(library(Seurat))

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
if (FALSE) { # debug
  opt = list(in_dir = "/home/lhou1/scRNAseq/DATA",  
             meta_file = "meta_info", sampleID_var = "patientID", cluster_var = "predicted.id",
             cluster = "EN-L2-3", group_var = "LesionType", reference_group = "AL", alternative_group = "CAL", count_file = "filter_counts",
             covars = 'Gender,DiedAge', covars_formula = "Gender+DiedAge", 
             minimum_percentage_cell = 0.05, out_dir =  "/home/mryals/scRNAseq/RESULTS/DEV/350kDE_REPORT", out_prefix = "New_350k")
}

if (is.na(opt$covars)) {
    covars = NULL 
} else {
    covars = unlist(strsplit(opt$covars, split = "[ ,;]+"))
}

reference_group = unlist(strsplit(opt$reference_group, split = "[ ,;]+"))

alternative_group = unlist(strsplit(opt$alternative_group, split = "[ ,;]+"))


# read the data

dat <- readRDS(file.path(opt$in_dir, "concat_not_scaled_raw_counts_only100geneCellFilter_witeMetadata_plusClinical_fixed_ordering.RDS"))
meta.data <- dat@meta.data
counts <- dat@assays$RNA@counts
dim(counts)
dim(meta.data)


meta.data$cell<-meta.data$cell_barcode
#meta.data <- meta.data #%>%
    #mutate(AL = ifelse(LesionType=='AL',1,0),
           #CAL = ifelse(LesionType=='CAL',1,0))


# delete missing clinical data cells
meta_info<-meta.data%>%dplyr:::filter(!is.na(patientID)) %>%
  dplyr::filter(  LesionType == "CAL" |  LesionType == "NWM" ) 
  

index = which(meta.data$cell %in% meta_info$cell_barcode)
filter_counts = counts[, index]

#rm(counts)

# Set up paths
prg.dir <- "/home/qwan/scRNAseq/PROGRAMS/DEV"


# Load the pipeline code
source(file.path(prg.dir, "pipeline_class_pseudoBulk_mode_ANCOVA.R"))
#rm(dat)
#rm(meta.data)

# Load the data
sce <- BiostatsSingleCell$new(count_data = filter_counts,
                              meta_data = meta_info,
                              sampleId_col = opt$sampleID_var,
                              cluster_col = opt$cluster_var,
                              treatment_col = opt$group_var)

# Make QC plots
#sce$make_QCplots(file.path(opt$out_dir, "UMI_diagnosis_QC_plot.pdf"))

# Filtering round 1
sce$apply_filter()

cluster = opt$cluster
prefix = opt$out_prefix

# set mode
sce$set_group_mode(cluster_of_interest = opt$cluster, ref_group = reference_group, alt_group = alternative_group)


# Filtering round 2
sce_qc <- sce$apply_filter_contrasts_R6()
#rm(sce)

# pull out filter info 

# filter<- sce$get_filter_info()
# filter1<-data_frame(data = filter) %>% 
#   group_by(name = names(data)) %>% 
#   unnest() %>%
#   mutate(i = row_number()) %>% 
#   spread(name, data)%>%dplyr:::select(-i)
# 
# write.csv(filter1, file.path(opt$out_dir, glue("{prefix}_filter_gene_info_{cluster}.csv")))


# ## Pull out normalized count 
# norm.count <- sce_qc$get_norm_counts() 
# ## Save the DE results to a file
# 
# norm.count <- cbind(rownames(norm.count), data.frame(norm.count, row.names=NULL))
# colnames(norm.count)[1]<-"V1"
# write_fst(norm.count,file.path(opt$out_dir,glue("{prefix}_norm_count_cell_{cluster}.fst") ) )


# ## Run t.test pipeline
t.test.res <- sce_qc$t_test_pipeline()
## Save the DE results to a file
write.csv(t.test.res, file.path(opt$out_dir, glue("{prefix}_t.test_cell_{cluster}.csv")), row.names = FALSE)
## Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "t-test")
## Save the volcano plot to a file
ggsave(file.path(opt$out_dir,glue("{prefix}_t.test_volcano_cell_{cluster}.png")))
# #
## Run u-test pipeline
u.test.res <- sce_qc$u_test_pipeline()
write.csv(u.test.res, file.path(opt$out_dir, glue("{prefix}_u.test_cell_{cluster}.csv")), row.names = FALSE)
## Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "u-test")
ggsave(file.path(opt$out_dir,glue("{prefix}_u.test_volcano_cell_{cluster}.png")))
# #
# ## Run pooled u-test
# #pooled.u.test.res <- sce_qc$pooled_u_test_pipeline()
# #write.csv(pooled.u.test.res, file.path(opt$out_dir, glue("{prefix}_pooled.u.test_cell_{cluster}.csv")), row.names = FALSE)
# ## Volcano plot
# #p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "Pooled u-test")
# #ggsave(file.path(opt$out_dir,glue("{prefix}_pooled.u.test_volcano_cell_{cluster}.png")))
#
# Run edgeR
edgeR.res <- sce_qc$edgeR_pipeline(covs = covars)#, covs_formula = opt$covars_formula)
write.csv(edgeR.res, file.path(opt$out_dir, glue("{prefix}_edgeR_cell_{cluster}.csv")), row.names = FALSE)
# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "edgeR")
ggsave(file.path(opt$out_dir,glue("{prefix}_edgeR_volcano_cell_{cluster}.png")))
#
# Run limma voom
limma.res <- sce_qc$limma_pipeline(covs = covars)
write.csv(limma.res, file.path(opt$out_dir, glue("{prefix}_limma_cell_{cluster}.csv")), row.names = FALSE)
# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "limma voom")
ggsave(file.path(opt$out_dir,glue("{prefix}_limma_volcano_cell_{cluster}.png")))
#
# Run DESeq2
# w/o shrink first
deseq2.res <- sce_qc$DESeq2_pipeline(covs = covars, shrink = FALSE)
write.csv(deseq2.res, file.path(opt$out_dir, glue("{prefix}_DESeq2_cell_{cluster}.csv")), row.names = FALSE)
# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "DESeq2 w/o shrink")
ggsave(file.path(opt$out_dir,glue("{prefix}_DESeq2_volcano_cell_{cluster}.png")))
# w shrink now
deseq2.shrink.res <- sce_qc$DESeq2_pipeline(covs = covars)
write.csv(deseq2.shrink.res, file.path(opt$out_dir, glue("{prefix}_DESeq2.shrink_cell_{cluster}.csv")), row.names = FALSE)
# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "DESeq2 with shrink")
ggsave(file.path(opt$out_dir,glue("{prefix}_DESeq2.shrink_volcano_cell_{cluster}.png")))


# Run glmmTMB
glmmTMB.res <- sce_qc$glmmTMB_pipeline(covs = covars, family = "nbinom2", detection_rate = FALSE)
write.csv(glmmTMB.res, file.path(opt$out_dir, glue("{prefix}_glmmTMB_cell_{cluster}_wo_cdr.csv")), row.names = FALSE)
# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "glmmTMB nbinom2")
ggsave(file.path(opt$out_dir,glue("{prefix}_glmmTMB_volcano_cell_{cluster}_wo_cdr.png")))
#
# Run MAST
MAST.res <- sce_qc$MAST_pipeline(covs = covars, detection_rate = TRUE)
write.csv(MAST.res, file.path(opt$out_dir, glue("{prefix}_MAST_cell_{cluster}_w_cdr.csv")), row.names = FALSE)
## Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "MAST")
ggsave(file.path(opt$out_dir,glue("{prefix}_MAST_volcano_cell_{cluster}_w_cdr.png")))
#

#
# Run limma cell level
limma.res <- sce_qc$limma_cell_level_pipeline(covs = covars)
write.csv(limma.res, file.path(opt$out_dir, glue("{prefix}_limma.cell.level_cell_{cluster}.csv")), row.names = FALSE)
# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "limma (cell level) voom")
ggsave(file.path(opt$out_dir,glue("{prefix}_limma.cell.level_volcano_cell_{cluster}.png")))
#
# Run ANCOVA
ancova.res <- sce_qc$ancova_pipeline(covs = covars)
write.csv(ancova.res, file.path(opt$out_dir, glue("{prefix}_ancova_cell_{cluster}.csv")), row.names=FALSE)
##Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "ancova")
ggsave(file.path(opt$out_dir, glue("{prefix}_ancova_volcano_cell_{cluster}.png")))


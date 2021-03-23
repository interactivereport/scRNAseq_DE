#.libPaths(c("/camhpc/pkg/R/3.5.1/centos6/lib64/R/library","/camhpc/home/jgagnon1/R/x86_64-pc-linux-gnu-library/3.5"))

suppressMessages(library(data.table))
suppressMessages(library(glue))
suppressMessages(library(scater))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(Matrix))
#suppressMessages(library(Seurat))

source("pipeline_class_042420.R")

option_list = list(
  make_option("--cluster", action = "store", default = "cell_type", type = "character",
              help = "Which cluster to use? [default: %default]"),
  make_option("--in_dir", action = "store", default = NA, type = "character",
              help = "Path for input file [required]"),
  make_option("--out_dir", action = "store", default = NA, type = "character",
              help = "Path for output file [required]"),
  make_option("--split", action = "store", default = 1, type = "integer",
              help = "partition to analyze [required]")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)

contrast = c("NWM","CAL")
my_rds = paste0(opt$in_dir,opt$cluster,opt$split,contrast[2],"vs",contrast[1],".rds")
no_filter <- BiostatsSingleCell$new(rds_file = my_rds,sampleId_col = "patientID",cluster_col = "predicted.id",treatment_col = "LesionType")
print("done loading.")

load(paste0(opt$in_dir,opt$cluster,contrast[2],"vs",contrast[1],"_meta_filtered.Rdata")) # my_libs
glmm_result = no_filter$glmmTMB_pipeline(input_lib_sizes=my_libs,cores=1,covs=c("Gender","DiedAge"),covs_formula = c("Gender","DiedAge"))

cluster = opt$cluster
split = opt$split
prefix = paste0(contrast[1],contrast[2])
write.csv(glmm_result, file.path(opt$out_dir, glue("{prefix}_glmmTMB_cell_{cluster}_{split}.csv")), row.names = FALSE)
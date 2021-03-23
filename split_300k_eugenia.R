source("pipeline_class_042420.R")

suppressMessages(library(data.table))
suppressMessages(library(glue))
suppressMessages(library(scater))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(Matrix))
suppressMessages(library(Seurat))

my_rds = "/home/kli3/proj/HTVC/MS/SingleCell/BatchCorrectionBenchMarking_tests/test6_all_Imperial_samples_pre_mRNA_ref/working_data/concat_not_scaled_raw_counts_only100geneCellFilter_witeMetadata_plusClinical_fixed_ordering.RDS"

dat <- readRDS(my_rds)
meta.data <- dat@meta.data
counts <- dat@assays$RNA@counts
meta.data$cell<-meta.data$cell_barcode
meta_info<-meta.data%>%dplyr:::filter(!is.na(patientID))
index = which(meta.data$cell %in% meta_info$cell_barcode)
filter_counts = counts[, index]
no_filter <- BiostatsSingleCell$new(count_data = filter_counts,
                              meta_data = meta_info,
                              sampleId_col = "patientID",
                              cluster_col = "predicted.id",
                              treatment_col = "LesionType")

filterR1 = no_filter$apply_filter_R6()

all_cell_types = unique(filterR1$pData()[,"predicted.id"])
all_cell_types = na.omit(all_cell_types)
print(all_cell_types)

ref_group = "NWM"
alt_group = c("CAL")

out_dir = paste0("/home/jgagnon1/jupyter-notebook-dir/",paste0(alt_group,collapse=""),"vs",ref_group,"_meta_filtered/")

for (cell_type in all_cell_types)
{
    print(cell_type)
    error_flag <- FALSE
    
    tryCatch({
        filterR1$set_group_mode(cell_type,ref_group,alt_group)
    }, error = function(err) {
         error_flag <<- TRUE    
        }
    )
    if (error_flag) next

    tryCatch({
        filterR2 = filterR1$apply_filter_contrasts_R6()
        }, error = function(err) {
         error_flag <<- TRUE
        }
    )
    if (error_flag) next
    
    ngenes = filterR2$get_genes()
    partition = split(seq(ngenes),cut_interval(seq(ngenes),16))
    my_pData = filterR2$pData()
    row.names(my_pData) = my_pData[,"cell"]
    for (split_num in seq(16))
    {
        obj = CreateSeuratObject(counts = filterR2$assayData()[partition[[split_num]],],project="J300k",meta.data=my_pData)
        saveRDS(obj,paste0(out_dir,cell_type,split_num,paste0(alt_group,collapse=""),"vs",ref_group,".rds"))
    }
    my_libs = filterR2$get_lib_sizes()
    save(my_libs,file=paste0(out_dir,cell_type,paste0(alt_group,collapse=""),"vs",ref_group,"_meta_filtered.Rdata"))
}
suppressMessages(library("EnsDb.Hsapiens.v86"))
suppressMessages(library(scater))
suppressMessages(library(tidyverse))
library(beachmat)

.libPaths("/camhpc/home/jgagnon1/R/x86_64-pc-linux-gnu-library/3.5/")

source("pipeline_class_042420.R")

opt = list(in_dir = "/home/jgagnon1/jupyter-notebook-dir", gene_info_file = "MS_Nature_gene_info.csv", 
          meta_file = "MS_Nature_meta_data.csv", sampleID_var = "sample", cluster_var = "cell_type",
          cluster = "Microglia", group_var = "diagnosis", reference_group = "Control", alternative_group = "MS", 
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

# load the data
sce <- BiostatsSingleCell$new(count_data = count_data,
                              meta_data = meta_info,
                              sampleId_col = opt$sampleID_var,
                              cluster_col = opt$cluster_var,
                              treatment_col = opt$group_var)
# Make QC plots
sce$make_QCplots(file.path(opt$out_dir, "UMI_diagnosis_QC_plot.pdf"))

# 1st round of filtering
sce$apply_filter(min.perc.cells.per.gene = 0.00) # 0% expression requirement

# set mode
sce$set_group_mode(cluster_of_interest = opt$cluster, ref_group = opt$reference_group, alt_group = opt$alternative_group)

# Filtering round 2
sce_qc <- sce$apply_filter_contrasts_R6()

cluster = opt$cluster
prefix = opt$out_prefix

# --------- t-test ------------#
# ## Run t.test pipeline
t.test.res <- sce_qc$t_test_pipeline()
write.csv(t.test.res, file.path(opt$out_dir, glue("{prefix}_t.test_cell_{cluster}.csv")), row.names = FALSE)

## Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "t-test")
ggsave(file.path(opt$out_dir,glue("{prefix}_t.test_volcano_cell_{cluster}.png"))) ## Save the volcano plot to a file


# --------- u-test -------------- #
## Run u-test pipeline
u.test.res <- sce_qc$u_test_pipeline()
write.csv(u.test.res, file.path(opt$out_dir, glue("{prefix}_u.test_cell_{cluster}.csv")), row.names = FALSE)

## Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "u-test")
ggsave(file.path(opt$out_dir,glue("{prefix}_u.test_volcano_cell_{cluster}.png")))

# ----------- edgeR -----------#
# Run edgeR
edgeR.res <- sce_qc$edgeR_pipeline(covs = covars)#, covs_formula = opt$covars_formula)
write.csv(edgeR.res, file.path(opt$out_dir, glue("{prefix}_edgeR_cell_{cluster}.csv")), row.names = FALSE)

# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "edgeR")
ggsave(file.path(opt$out_dir,glue("{prefix}_edgeR_volcano_cell_{cluster}.png")))

# ------------- limma-voom -------- #

# Run limma voom
limma.res <- sce_qc$limma_pipeline(covs = covars)
write.csv(limma.res, file.path(opt$out_dir, glue("{prefix}_limma_cell_{cluster}.csv")), row.names = FALSE)

# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "limma voom")
ggsave(file.path(opt$out_dir,glue("{prefix}_limma_volcano_cell_{cluster}.png")))

# ------- DEseq2 -------------- #
# with shrink 
deseq2.shrink.res <- sce_qc$DESeq2_pipeline(covs = covars)
write.csv(deseq2.shrink.res, file.path(opt$out_dir, glue("{prefix}_DESeq2.shrink_cell_{cluster}.csv")), row.names = FALSE)

# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "DESeq2 with shrink")
ggsave(file.path(opt$out_dir,glue("{prefix}_DESeq2.shrink_volcano_cell_{cluster}.png")))

# --------- glmmtmb ------------- #
# Run glmmTMB
glmmTMB.res <- sce_qc$glmmTMB_pipeline(covs = covars, family = "nbinom2", detection_rate = FALSE)
write.csv(glmmTMB.res, file.path(opt$out_dir, glue("{prefix}_glmmTMB_cell_{cluster}_wo_cdr.csv")), row.names = FALSE)

# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "glmmTMB nbinom2")
ggsave(file.path(opt$out_dir,glue("{prefix}_glmmTMB_volcano_cell_{cluster}_wo_cdr.png")))

# -------- MAST ----------------- #
# Run MAST
MAST.res <- sce_qc$MAST_pipeline(covs = covars, detection_rate = TRUE)
write.csv(MAST.res, file.path(opt$out_dir, glue("{prefix}_MAST_cell_{cluster}_w_cdr.csv")), row.names = FALSE)

## Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "MAST")
ggsave(file.path(opt$out_dir,glue("{prefix}_MAST_volcano_cell_{cluster}_w_cdr.png")))

#---------------------
# Run limma cell level
# -------------------
limma.cell.res <- sce_qc$limma_cell_level_pipeline(covs = covars)
write.csv(limma.cell.res, file.path(opt$out_dir, glue("{prefix}_limma.cell.level_cell_{cluster}.csv")), row.names = FALSE)

# Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "limma (cell level) voom")
ggsave(file.path(opt$out_dir,glue("{prefix}_limma.cell.level_volcano_cell_{cluster}.png")))

#-----------------
# Run ANCOVA
# ---------------
ancova.res <- sce_qc$ancova_pipeline(covs = covars)
write.csv(ancova.res, file.path(opt$out_dir, glue("{prefix}_ancova_cell_{cluster}.csv")), row.names=FALSE)

##Volcano plot
p <- sce_qc$volcanoPlot(FDR_threshold = 0.05, FC_threshold = 2, title = "ancova")
ggsave(file.path(opt$out_dir, glue("{prefix}_ancova_volcano_cell_{cluster}.png")))

# -------- Violin plots ----------------- #
# t-test violin
p <- sce_qc$violinPlot(gene.name="SLC1A3",de.method="t.test",cellinfo="cell")
ggsave(file.path(opt$out_dir,glue("{prefix}_t.test_violin_cell_{cluster}_SLC1A3.png")))


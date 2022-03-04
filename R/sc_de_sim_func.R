#' @title sc_de_sim_func.R
#' @description Generate Single Cell DE Simulation Results
#' @details This function runs a single cell DE simulation to obtain DE method performance for DE methods
#' @param in_dir The directory where the input files are located in 10X cellranger v2 format
#' @param gene_info_file The filename of the gene info file, csv format
#' @param gene_name_column The name of the column containing unique gene identifiers in the gene info file
#' @param meta_file The filename of the cell metadata file, csv format
#' @param sample_column The name of the sample column in the metadata file
#' @param cell_type_column The name of the cell type in the metadata file
#' @param cell_barcode_column The name of the cell barcode indices column in the metadata file
#' @param contrast_column The name of the DE contrast column in the metadata file, needs 2 levels, needs to map to samples
#' @param contrast_ref_group The level of the contrast to use as a reference
#' @param contrast_alt_group The level of the contrast to use as an alternative group
#' @param count_file The name of the count file (matrix market format)
#' @param cluster The name of the cell type to use for the contrast
#' @param covars The covariates in the metadata file
#' @param simulation_type Simulation type to run.  Values allowed are 'de' and 'null'.  Default is 'de'.
#' @param fc_estimate If running a when `simulation_type` is set to 'de', set `fc_estimate` to TRUE to run.  Default is FALSE.
#' @param simulation_mode Change the way the simulation is performed.  Values allowed are 1,2, and 3.  Default is 3.
#' @param fc The simulated fold change.  When `simulation_type` set to 'null' then simulations `fc` is set to 1.
#' @param pDE The percentage of genes that will be DE in the simulation
#' @param nSim.subjects The number of subjects to simulate in each contrast group
#' @param sim_seed_mode Mode to use for seed.  'nsim_mode' will generate seeds based on the `n_simulation`.  'seed_mode' will run the seeds specified by the user by `simulation_seeds`.
#' Default value is 'nsim_mode'.
#' @param n_simulation The number of simulations to perform.  Default is 1.
#' @param simulation_seeds A vector of simulation seeds to run.
#' @param main_seed Set a seed for reproducible simulations.  This seed is used in both modes of `sim_seed_mode`.
#' @param percentage.cell Percentage of cells to downsample for use in the simulation.  Useful for cell types with large numbers of cells.  Default is 1, i.e. no downsampling.
#' @param minimum.cells.per.gene.type Contrast-level filtering.  Values allowed are 'or' and 'and'.  Default is 'or'.  See public method apply_filter_contrasts_R6 in BiostatsSingleCell for details.
#' @param minimum.cells.per.subject Contrast-level filtering.  Specifies minimum number of cells per subject to retain.  Default is 5.  See public method apply_filter_contrasts_R6 in BiostatsSingleCell for details.
#' @param sim_cores Number of cores to use to run the simulations.  Default is 5.
#' @param de_methods DE methods to run from. Allowed values include 'ancova', 't_test', 'u_test', 'edgeR', 'limma', 'DESeq2.shrink', 'DESeq2', 'DEseq2.nofilt', 'nebula_LN', 'nebula_HL', 'glmmTMB', 'MAST'
#' Please note when including 'MAST' in de_methods, you will only be able to run one simulation at a time in either mode.
#' Default value is all methods excluding MAST.
#' @details The `simulation_mode` argument runs with three different modes
#' \itemize{
#' \item Mode 1 : Simulate mean/disp. in ref and alt groups based on mean/disp. from real data ref and alt groups, respectively.
#' \item Mode 2 : Simulate ref mean/disp. from real ref data; derive alt mean estimate using real ref mean multiplied by simulation FC;
#' derive alt disp. estimate either from real alt data for abundant genes or from median disp. for lowly-expressed genes.
#' \item Mode 3 : Simulate mean/disp. from real ref data; derive alt mean estimate using real ref mean multiplied by simulation FC;
#' predict alt disp. by testing simulation alt. mean in general additive model (GAM) after training the GAM on disp. ~ mean from real data
#' }
#' @return A dataframe of relevant simulation data (FC bias, null alpha, or DE performance parameters)
#' @import glmmTMB
#' @import Matrix
#' @import SingleCellExperiment
#' @import dplyr
#' @import purrr
#' @import fitdistrplus
#' @import Hmisc
#' @import edgeR
#' @import gridExtra
#' @import SimSeq
#' @import fdrtool
#' @import ROCR
#' @import tidyr
#' @import PRROC
#' @import data.table
#' @import mgcv
#' @import parallel
#' @import MASS
#' @export

sc_de_sim_func <- function(in_dir,gene_info_file,gene_name_column,meta_file,sample_column,
                           cell_type_column,cell_barcode_column,contrast_column,
                           contrast_ref_group,contrast_alt_group,count_file,
                           cluster,covars,simulation_type = "de",fc_estimate = FALSE,
                           simulation_mode = 3,fc,pDE,nSim.subjects,n_simulation=1,
                           main_seed,percentage.cell=1,minimum.cells.per.gene.type = "or",
                           minimum.cells.per.subject = 5,sim_cores = 5,
                           de_methods = c('ancova', 't_test', 'u_test', 'edgeR', 'limma', 'DESeq2.shrink',
                                      'DESeq2', 'DEseq2.nofilt', 'nebula_LN', 'nebula_HL', 'glmmTMB'),
                           sim_seed_mode = "nsim_mode", simulation_seeds
                       ){

  #Take care of bindings for global variables
  Var2 = NULL
  Freq = NULL
  my_disp = NULL
  ID = NULL
  tFC = NULL
  log2FC = NULL
  sFC = NULL


  # user must enter a cell ranger v2 format data
  count_data <- readMM(file.path(in_dir,count_file))
  meta_info=read.csv(file.path(in_dir, meta_file), stringsAsFactors=FALSE)
  gene_info=read.csv(file.path(in_dir, gene_info_file), stringsAsFactors=FALSE)

  #Data format checks - can we think of any others
    #Is your gene identifier colum in your gene info file
    if(!gene_name_column %in% colnames(gene_info)){
      stop("Gene name column not present in gene info file.")
    }
    #Is your cell barcode column in your metadata?
    if(!cell_barcode_column %in% colnames(meta_info)){
      stop("Cell barcode column not present in the metadata file.")
    }

    #Pipeline requires barcoded column named 'cell'
    meta_info$cell = meta_info[,cell_barcode_column]

    #Is your contrast column in your metadata?
    if(!contrast_column %in% colnames(meta_info)){
      stop("Contrast column is not present in the metadata file.")
    }
    #Is your cell type column in your metadata?
    if(!cell_type_column %in% colnames(meta_info)){
      stop("Cell type column is not present in the metadata file.")
    }
    #Is your sample column in your metadata?
    if(!sample_column %in% colnames(meta_info)){
      stop("Sample column does not exist in your metadata file")
    }
    #Are all of your covariate columns in your metadata?
    if(!all(covars %in% colnames(meta_info))){
      stop("Covariate columns do not exist in your metadata file.")
    }
    #Are the dimensions of your gene info the same as your counts?
    if(nrow(gene_info) != nrow(count_data)){
      stop("The number of rows in your gene info does not match number of count rows.")
    }
    #Are the dimensions of your metadata the same as your counts?
    if(nrow(meta_info) != ncol(count_data)){
      stop("The number of rows in your metadata does not match number of count columns.")
    }

    #Pipeline formatting checks
    #Does your contrast column have two levels?
    if(length(levels(as.factor(meta_info[,contrast_column]))) != 2){
      print(table(meta_info[,contrast_column], useNA = "always"))
      print(length(levels(as.factor(addNA(meta_info[,contrast_column])))))
      stop("Contrast column does not contain exactly two levels.")
    }
    #Is your reference group in the contrast column?
    if(!contrast_ref_group %in% unique(meta_info[,contrast_column])){
      stop("Contrast reference group not found in contrast column.")
    }
    #Is your alternative group in the contrast column?
    if(!contrast_alt_group %in% unique(meta_info[,contrast_column])){
      stop("Contrast alternative group not found in contrast column.")
    }
    #Is your cluster in your cell type column?
    if(!cluster %in% unique(meta_info[,cell_type_column])){
      stop("Cluster is not present in cell type column.")
    }

  #Do a check of simulation mode and warnings
  if(fc_estimate){
    if(simulation_type == "null"){
      warning("fc_estimate TRUE is only calculated for the simulation type 'de'.")
    }
  }
  #Set fc to 1 if in null mode
  if(simulation_type == "null"){
    fc <- 1
    warning("simulation type is null, overriding argument 'fc' and setting 'fc' to 1.")
  }

  #Do a check of de methods to handle MAST
  if("MAST" %in% de_methods){
    #MAST requires 1 simulation only
    if(sim_seed_mode == "seed_mode"){
      if(length(unique(simulation_seeds)) != 1){
          stop("simulation_seeds must be set to a single seed value when 'MAST' is selected as a de_method.")
        }
      }
    if(sim_seed_mode == "nsim_mode"){
      if(n_simulation != 1){
          stop("n_simulation must be set to 1 when 'MAST' is selected as a de_method.")
        }
      }
    warning("You are performing a MAST simulation.  MAST is slower than other methods.")
  }

  colnames(count_data) <- meta_info[,cell_barcode_column]
  rownames(count_data) <- gene_info[,gene_name_column]

  # assume equal subjects for Control/MS
  nControlSim = nSim.subjects # number of Controls to simulate
  nMSSim = nSim.subjects # number of MS to simulate

  # determine Control and MS samples
  my_samples = unique(meta_info[,sample_column])
  my_cell_types = unique(meta_info[,cell_type_column])

  samples_table = data.frame(table(meta_info[,sample_column],meta_info[,contrast_column]))

  #Control/MS naming convention for ref/alt
  ControlSamples = as.character((samples_table %>% dplyr::filter(Var2 == contrast_ref_group,Freq > 0))$"Var1")
  MSSamples = as.character((samples_table %>% dplyr::filter(Var2 != contrast_ref_group,Freq > 0))$"Var1")

  stopifnot(nControlSim == nMSSim)
  stopifnot(nControlSim <= length(ControlSamples))
  stopifnot(nMSSim <= length(MSSamples))

  # set mode
  mode = simulation_mode

  # loop over cell types
  # The piece of code below this is not a loop until the lapply statement
  acelltype = cluster

  idx = which(meta_info[, cell_type_column] == acelltype)

  count_subset = count_data[,idx]
  meta_subset = meta_info[idx,]

  # filter out low expressed genes
  count_sums = apply(count_subset,1,function(x) sum(x > 0)) #LP

  rows_idx = which(count_sums == 0)
  count_subset = count_subset[-rows_idx,]  ## remove genes with all 0 counts across cells
  gene_info_subset =  gene_info[-rows_idx,]

  # setup DE genes and fcs
  ngenes = nrow(count_subset)
  set.seed(main_seed)
  DE_idx = sample(seq(ngenes),floor(ngenes*pDE),replace=FALSE)

  DE_list <- gene_info_subset[DE_idx, gene_name_column]
  NDE_idx = setdiff(seq(ngenes),DE_idx)
  fc_vec = rep(1,ngenes)

  set.seed(main_seed*10)
  fc_vec[DE_idx] = sample(c(fc,1/fc),length(DE_idx),replace = TRUE)

  # add data set to save the gene and true FC
  if(fc_estimate){
    fcdat2 <- data.frame("ID"=rep(NA,ngenes),"tFC"=rep(NA,ngenes))
    fcdat2$ID <- gene_info_subset[,gene_name_column]
    fcdat2$tFC <- fc_vec
  }

  # estimate mean and disp for each sample
  my_mu_disp = lapply (my_samples, function(asample)
  {
    print(asample)
    #sample_idx = meta_subset$sample == asample
    sample_idx <- which(meta_subset[,sample_column] == asample)
    # if there is 1 or few cells return NULL
    if (sum(sample_idx) <= 1) return(list("avCPM"=NULL,"disp"=NULL,"libs"=NULL,"ncell"=0,"meta"=NULL))
    counts_sample = count_subset[,sample_idx]
    meta_sample = meta_subset[sample_idx,]

    # find rows with at least 1 cell with UMI > 1
    counts_sum = apply(counts_sample,1,function(x) sum(x > 0)) #LP
    good_counts = which(counts_sum > 0)

    # disp is 0 by default
    tag_disp = rep(0,nrow(counts_sample))

    # fit model only to good rows
    y <- DGEList(counts_sample[good_counts,])
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y) # !! might be affected by large number of zeros in sample level data !!
    y <- estimateTagwiseDisp(y) # dispersion calculated at sample level, so capbatch, sex, age is fixed

    # good rows get a disp value, bad rows get disp = 0
    tag_disp[good_counts] = y$tagwise.dispersion

    cpm_Jake = sweep(counts_sample,2,y$samples$lib.size,"/")*10^6  ## library size normalization

    temp_avCPM = rep(0,nrow(counts_sample))
    temp_avCPM[good_counts] = apply(cpm_Jake[good_counts,],1,mean)

    # good rows are in good_idx
    list("avCPM" = temp_avCPM, "disp"=tag_disp,"libs"=y$samples$lib.size,"ncell"=ncol(counts_sample),"meta"=meta_sample,"good_idx"=good_counts)
  }
  )
  names(my_mu_disp) = my_samples

  #Simulation function - calls R6 pipeline instance
  run_simulation <- function(){

    my_Csamples = sample(ControlSamples, nControlSim,replace=FALSE)
    my_MSsamples = sample(MSSamples,nMSSim,replace=FALSE)

    mode = simulation_mode

    print("generating Ccells...")
    Ccells_list = lapply (my_Csamples, function(CSample){
      mu_disp = my_mu_disp[[CSample]]
      n_cells = mu_disp$ncell
      if (n_cells == 0) return(NULL)

      # my default simulated data is 0
      my_results = matrix(0,ngenes,n_cells)

      # use only good rows for the mu and theta calculation
      mu_vec = as.vector(outer(mu_disp$avCPM[mu_disp$good_idx]/10^6, mu_disp$libs)) #libs ok
      theta_vec = as.vector(outer(1/mu_disp$disp[mu_disp$good_idx],rep(1,n_cells)))

      # good rows are simulated from negbin, badrows are set to 0
      my_results[mu_disp$good_idx,]=matrix(rnegbin(n_cells*length(mu_disp$good_idx),mu=mu_vec,theta = theta_vec),nrow=length(mu_disp$good_idx),ncol=n_cells)
      my_results
    })
    Ccells = do.call("cbind",Ccells_list)

    print("generating MS cells...")
    MScell_list = lapply(seq(nMSSim),function(MS_idx){
      mu_dispMS = my_mu_disp[[my_MSsamples[MS_idx]]]
      mu_dispC = my_mu_disp[[my_Csamples[MS_idx]]]
      n_cells = mu_dispMS$ncell
      if (n_cells == 0) return(NULL)
      if (mode >= 2) if(mu_dispC$ncell == 0) return(NULL)

      # my default simulated data is 0
      my_results = matrix(0,ngenes,n_cells)

      # use good rows to calculate mean
      if (mode == 1) mu_vec = as.vector(outer(mu_dispMS$avCPM[mu_dispMS$good_idx]/10^6, mu_dispMS$libs))
      if (mode >= 2) mu_vec = as.vector(outer(mu_dispC$avCPM[mu_dispC$good_idx]/10^6*fc_vec[mu_dispC$good_idx], mu_dispMS$libs))

      # use good rows to calculate theta
      if (mode == 1)  theta_vec = as.vector(outer(1/mu_dispMS$disp[mu_dispMS$good_idx],rep(1,n_cells)))

      if (mode == 2)
      {
        # median dispersion
        med_disp = median(mu_dispMS$disp[mu_dispMS$good_idx])

        # default value is 0
        theta_vec = rep(0,ngenes)

        #genes which are expressed in both control and MS
        both_good = intersect(mu_dispC$good_idx,mu_dispMS$good_idx)
        theta_vec[both_good]=mu_dispMS$disp[both_good]

        # genes expressed in control but not in MS
        MS_bad_idx = setdiff(seq(ngenes),mu_dispMS$good_idx)
        control_good = intersect(mu_dispC$good_idx,MS_bad_idx)
        theta_vec[control_good] = med_disp # MS expression is not good, so use imputed value

        theta_vec = theta_vec[mu_dispC$good_idx]
      }
      if (mode == 3)
      {
        # use good rows to estimate mean/disp relationship
        temp_mean = mu_dispMS$avCPM[mu_dispMS$good_idx]/10^6
        temp_data = data.frame("my_mean" = temp_mean,"my_disp" = mu_dispMS$disp[mu_dispMS$good_idx])

        # filter out very low dispersions
        filter_data = temp_data %>% dplyr::filter(my_disp > min(temp_data$my_disp))

        # fit mean/disp relationship
        my_gam=gam(my_disp~s(my_mean,bs="cs"),data=filter_data)

        my_predict = pmax(predict(my_gam,newdata=data.frame("my_mean"=mu_dispC$avCPM[mu_dispC$good_idx]/10^6*fc_vec[mu_dispC$good_idx])),min(mu_dispMS$disp[mu_dispMS$good_idx]))
        theta_vec = as.vector(outer(1/my_predict,rep(1,n_cells)))
      }
      if (mode == 1) my_results[mu_dispMS$good_idx,]=matrix(rnegbin(n_cells*length(mu_dispMS$good_idx),mu=mu_vec,theta=theta_vec),nrow=length(mu_dispMS$good_idx),ncol=n_cells)
      if (mode >= 2)
      {
        my_results[mu_dispC$good_idx,]=matrix(rnegbin(n_cells*length(mu_dispC$good_idx),mu=mu_vec,theta=theta_vec),nrow=length(mu_dispC$good_idx),ncol=n_cells)

        # if row is DE and bad expression in control but good expression in MS:
        temp_idx = intersect(intersect(DE_idx, setdiff(seq(ngenes),mu_dispC$good_idx)),mu_dispMS$good_idx)
        temp_mu = as.vector(outer(mu_dispMS$avCPM[temp_idx]/10^6, mu_dispMS$libs))
        temp_theta = as.vector(outer(1/mu_dispMS$disp[temp_idx],rep(1,n_cells)))
        my_results[temp_idx,] = matrix(rnegbin(n_cells*length(temp_idx),mu=temp_mu,theta=temp_theta),nrow=length(temp_idx),ncol=n_cells)


        print(length(intersect(intersect(DE_idx, setdiff(seq(ngenes),mu_dispC$good_idx)),setdiff(seq(ngenes),mu_dispMS$good_idx))))

      }
      my_results
    })
    # bind simulated data
    MScells = do.call("cbind",MScell_list)
    cluster_sim = cbind(Ccells,MScells)
    cluster_sim<-as(cluster_sim, "dgTMatrix")

    # Control meta data
    Cmeta = lapply(my_Csamples,function(CSample) my_mu_disp[[CSample]]$meta) # meta is already null if ncells = 0

    # meta data for MS samples
    if (mode == 1) MSmeta = lapply(my_MSsamples,function(MSSample) my_mu_disp[[MSSample]]$meta)
    if (mode >= 2)
    {
      # use Control meta information, map Control sample labels to MS sample labels
      MSmeta = lapply(my_MSsamples,function(MSSample) my_mu_disp[[MSSample]]$meta)

      # replace alt meta data age,sex and capbatch by Control meta data
      #Confounding covars are not simulated, map control metadata to alt metadata
      for(i in 1:length(MSmeta)){
        MSmeta_select <- MSmeta[[i]]
        Cmeta_select <- Cmeta[[i]]
        for(j in covars){
          MSmeta_select[,j] <- rep( unique(Cmeta_select[,j]), length(MSmeta_select[,sample_column]) ) #can the sample column be used instead of 'age'
        }
        MSmeta[[i]] <- MSmeta_select
      }

      cluster_meta = rbind(do.call("rbind",Cmeta),do.call("rbind",MSmeta))

    }

    # all meta data
    cluster_meta = rbind(do.call("rbind",Cmeta),do.call("rbind",MSmeta))

    # add gene name as rownames and cell name as colnames in count data
    rownames(cluster_sim) <- rownames(count_subset) #gene_info$index
    #colnames(cluster_sim) <- cluster_meta$cell
    colnames(cluster_sim) <- cluster_meta[,cell_barcode_column]

    # downsample cell in each subject
    cluster_meta <- cluster_meta %>% dplyr::group_by(!!sym(sample_column)) %>%
      dplyr::sample_frac(percentage.cell) %>% dplyr::ungroup() %>% as.data.frame()

    # subset count data column barcode names to the downsampled cell barcodes
    cluster_sim<-cluster_sim[,cluster_meta[,cell_barcode_column]]

    sce <- BiostatsSingleCell$new(count_data = cluster_sim,
                                  meta_data = cluster_meta,
                                  sampleId_col = sample_column,
                                  cluster_col = cell_type_column,
                                  treatment_col = contrast_column)

    # Filtering round 1
    sce$apply_filter(min.perc.cells.per.gene = 0.00)

    cluster = acelltype

    # set mode
    sce$set_group_mode(cluster_of_interest = cluster, ref_group = contrast_ref_group, alt_group = contrast_alt_group)

    # Filtering round 2
    sce_qc <- sce$apply_filter_contrasts_R6(min.cells.per.gene.type=minimum.cells.per.gene.type,
                                            min.perc.cells.per.gene=0.1,
                                            min.cells.per.subj=minimum.cells.per.subject)

    #Internal performance function
    get_performance <- function(res, id = "ID", pval = "Pvalue", padj = 'FDR',
                                de.genes = "", fc = "log2FC", fc.threshold = 1.1, use.FC = FALSE){
      # remove genes if padj == NA
      res <- res[!is.na(res[, padj]), ]
      res$DE.known <- ifelse(res[, id] %in% de.genes, TRUE, FALSE)
      if (use.FC) {
        log2fc <- log2(fc.threshold)
        res$DE.test <- ifelse(res[, padj] < 0.05 & abs(res[, fc]) > log2fc, TRUE, FALSE)
      } else {
        res$DE.test <- ifelse(res[, padj] < 0.05, TRUE, FALSE)
      }
      power.res <- sum(res$DE.known & res$DE.test)/sum(res$DE.known)
      fdr.res <- sum(!res$DE.known & res$DE.test)/sum(res$DE.test)

      # Convert p values to z score for auc calculation
      np <- length(which(res[, pval] < 0 | res[, pval] > 1))
      if (np > 0) {
        warning(np, " invalid p-values set to NA")
        res[, pval][res[, pval] < 0 | res[, pval] > 1] <- NA
      }

      res$z.score <- -qnorm((res[, pval]/2), FALSE)
      res$z.score[!is.finite(res$z.score)] <- NA

      pred <- try(prediction(predictions = res$z.score, res$DE.known), silent = TRUE)
      if (inherits(pred, "try-error")) {
        auc.res <- NA
        ppauc.res <-NA
      } else {
        auc.res <- performance(pred, measure = "auc")@y.values[[1]]
        pred.df <- data.frame("predictions"=pred@predictions[[1]], "labels"=pred@labels[[1]], stringsAsFactors = F)
        fg2 <- pred.df$predictions[pred.df$labels == TRUE]
        bg2 <- pred.df$predictions[pred.df$labels == FALSE]

        pr <- pr.curve(scores.class0 = fg2, scores.class1 = bg2, curve = T)
        prauc.res <- pr[[2]]

      }
      return(list(power = power.res, fdr = fdr.res, auc = auc.res,prauc=prauc.res))
    }

    # run each pipeline in de_methods

    res_list <- list()

    time_list <- list()

    if("ancova" %in% de_methods){

      time_ancova <- system.time({ancova_results = sce_qc$ancova_pipeline(covs = covars)})
      ancova_results <- ancova_results[!is.na(ancova_results$Pvalue), ]

      res_list$ancova_results <- ancova_results
      time_list$time_ancova <- time_ancova

    }

    if("t_test" %in% de_methods){

      time_t.test <- system.time({t_test_results = sce_qc$t_test_pipeline()})
      t_test_results <- t_test_results[!is.na(t_test_results$Pvalue), ]

      res_list$t_test_results <- t_test_results
      time_list$time_t.test <- time_t.test

    }

    if("u_test" %in% de_methods){

      time_u.test <- system.time({u_test_results =sce_qc$u_test_pipeline()})
      u_test_results <- u_test_results[!is.na(u_test_results$Pvalue), ]

      res_list$u_test_results <- u_test_results
      time_list$time_u.test <- time_u.test

    }

    if("edgeR" %in% de_methods){

      time_edgeR <- system.time({edgeR_results= sce_qc$edgeR_pipeline(covs = covars)})
      edgeR_results <- edgeR_results[!is.na(edgeR_results$Pvalue), ]

      res_list$edgeR_results <- edgeR_results
      time_list$time_edgeR <- time_edgeR

    }

    if("limma" %in% de_methods){

      time_limma <- system.time({limma_results = sce_qc$limma_pipeline(covs = covars)})
      limma_results <- limma_results[!is.na(limma_results$Pvalue), ]

      res_list$limma_results <- limma_results
      time_list$time_limma <- time_limma

    }

    if("DESeq2.shrink" %in% de_methods){

      time_DESeq2.shrink <- system.time({DESeq2.shrink_results = sce_qc$DESeq2_pipeline(covs = covars)})
      DESeq2.shrink_results <- DESeq2.shrink_results[!is.na(DESeq2.shrink_results$Pvalue), ]

      res_list$DESeq2.shrink_results <- DESeq2.shrink_results
      time_list$time_DESeq2.shrink <- time_DESeq2.shrink

    }

    if("DESeq2" %in% de_methods){

      time_DESeq2 <- system.time({DESeq2_results = sce_qc$DESeq2_pipeline(covs = covars,shrink = FALSE)})
      DESeq2_results <- DESeq2_results[!is.na(DESeq2_results$Pvalue), ]

      res_list$DESeq2_results <- DESeq2_results
      time_list$time_DESeq2 <- time_DESeq2

    }

    if("DEseq2.nofilt" %in% de_methods){

      time_DESeq2.nofilt <- system.time({DESeq2.nofilt_results = sce_qc$DESeq2_pipeline(covs = covars,shrink = FALSE, independent_filtering = FALSE)})
      DESeq2.nofilt_results <- DESeq2.nofilt_results[!is.na(DESeq2.nofilt_results$Pvalue), ]

      res_list$DESeq2.nofilt_results <- DESeq2.nofilt_results
      time_list$time_DESeq2.nofilt <- time_DESeq2.nofilt

    }

    if("nebula_LN" %in% de_methods){

      time_nebula.LN <- system.time({nebula.LN_results = sce_qc$nebula_pipeline(covs = covars,method="LN")$res.tab})
      nebula.LN_results <- nebula.LN_results[!is.na(nebula.LN_results$Pvalue),]

      res_list$nebula.LN_results <- nebula.LN_results
      time_list$time_nebula.LN <- time_nebula.LN

    }

    if("nebula_HL" %in% de_methods){

      time_nebula.HL <- system.time({nebula.HL_results = sce_qc$nebula_pipeline(covs = covars,method="HL")$res.tab})
      nebula.HL_results <- nebula.HL_results[!is.na(nebula.HL_results$Pvalue),]

      res_list$nebula.HL_results <- nebula.HL_results
      time_list$time_nebula.HL <- time_nebula.HL

    }

    if("glmmTMB" %in% de_methods){

      time_glmmTMB <- system.time({glmmTMB_results = sce_qc$glmmTMB_pipeline(covs = covars, family = "nbinom2",cores = 4 )})
      glmmTMB_results <- glmmTMB_results[!is.na(glmmTMB_results$Pvalue), ]
      nTotal <- nrow(glmmTMB_results)
      nError <- sum(glmmTMB_results$Error, na.rm = T)
      glmmTMB_percError <- paste0(round(nError/nTotal*100, 2), "%")

      res_list$glmmTMB_results <- glmmTMB_results
      time_list$time_glmmTMB <- time_glmmTMB

    }

    if("MAST" %in% de_methods){

      time_MAST <- system.time({MAST_results = sce_qc$MAST_pipeline(covs = covars,detection_rate = TRUE)})
      MAST_results <- MAST_results[!is.na(MAST_results$Pvalue), ]

      res_list$MAST_results <- MAST_results
      time_list$time_MAST <- time_MAST

    }

    if(simulation_type == "null"){

      #Do alpha calc
      alpha_vec <- unlist( lapply(res_list, function(x){ sum(x$Pvalue < 0.05) / nrow(x)}) )
      names(alpha_vec) <- gsub("results","alpha",names(alpha_vec))

      #Make res
      null_res <- c(n.gene.r2 = length(sce_qc$get_genes()),
                    n.cell = length(sce_qc$get_cells()),
                    alpha_vec
                    )

      return(null_res)

    }

    if(simulation_type == "de"){

      #Do performance
      perf_list <- lapply(res_list, get_performance, de.genes=DE_list)

      power_vec <- unlist(lapply(perf_list, "[[", "power"))
      names(power_vec) <- gsub("results","power",names(power_vec))

      fdr_vec <- unlist(lapply(perf_list, "[[", "fdr"))
      names(fdr_vec) <- gsub("results","fdr",names(fdr_vec))

      auc_vec <- unlist(lapply(perf_list, "[[", "auc"))
      names(auc_vec) <- gsub("results","auc",names(auc_vec))

      prauc_vec <- unlist(lapply(perf_list, "[[", "prauc"))
      names(prauc_vec) <- gsub("results","prauc",names(prauc_vec))

      time_vec <- unlist(lapply(time_list, "[[", "elapsed"))

      #Make res
      de_res <- c(
        n.gene.r2 = length(sce_qc$get_genes()),
        n.de.r2 = length( intersect(DE_list, sce_qc$get_genes() ) ),
        n.cell = length(sce_qc$get_cells()),
        power_vec,
        fdr_vec,
        auc_vec,
        prauc_vec,
        time_vec
      )

      if(fc_estimate){
        #Create the fc object
        res_list_select <- lapply(res_list, "[", c("ID","log2FC","Pvalue","FDR"))

        fcdat3_names <- gsub("\\_results","",names(res_list_select))

        names(res_list_select) <- fcdat3_names

        res_select_df <- res_list_select %>%
          purrr::imap(function(x, y) x %>% rename_with(~paste(., y, sep = '_'), -ID)) %>%
          purrr::reduce(full_join, by = 'ID') #retain full join list of gene IDs from all methods

        fcdat3 <- fcdat2 %>% left_join(res_select_df) #left_join here to join to fcdat2

        fcdat4 <- fcdat3 %>%
            dplyr::filter(ID %in% sce_qc$get_genes()) %>%
            stats::na.omit() %>%
            dplyr::select(ID, tFC, contains("log2FC")) %>%
            dplyr::filter(tFC != 1) %>%
            dplyr::rename_at(vars(matches("log2FC")), funs(gsub("log2FC\\_","",.))) %>%
            tidyr::pivot_longer(cols = c(-ID, -tFC), names_to="method", values_to="log2FC", values_drop_na = TRUE) %>%
            dplyr::mutate(sFC =2^log2FC,BIAS=abs(sFC-tFC)) #%>%
            #dplyr::filter(FC<10)

          return(list(de_res, fcdat4))
        }

      return(de_res)
  }

  }

  #Run simulation(s)
  #Check simulation mode
  if(!sim_seed_mode %in% c("nsim_mode", "seed_mode")){
    stop("sim_seed_mode must be 'nsim_mode' or 'seed_mode'.")
  }

  if(sim_seed_mode == "nsim_mode"){
  # Set seed
  set.seed(main_seed)
  seeds <- sample.int(1000000, n_simulation)
  message(paste0("Seeds for ", n_simulation, " simulations: ", paste(seeds, collapse=",")))
  }

  if(sim_seed_mode == "seed_mode"){
    seeds <- simulation_seeds
  }

  sim.res <- mclapply(seeds, function(seed){
    set.seed(seed)
    run_simulation()
  }, mc.cores = sim_cores) #set sim_cores to 1 is lapply

  return(sim.res)

}

#!/bin/bash

cell_names="Astrocytes B_cells Endo_cells EN-L2-3 EN-L4 EN-L5-6 EN-MIX EN-PYR Glia-MIX IN-PVALB IN-SST IN-SV2C IN-VIP Microglia OL OPC Phagocytes Stromal_cells T_cells"

#cell_names="EN-L2-3 EN-L4 EN-L5-6"

for my_cell_type in $cell_names
do
	echo "source /etc/profile.d/modules_bash.sh; module load R/3.6.1;Rscript explore_downsample.R --in_dir /home/jgagnon1/jupyter-notebook-dir --gene_info_file MS_Nature_gene_info.csv --count_file MS_Nature_UMI_subset.mtx --meta_file MS_Nature_meta_data.csv --sampleID_var sample --cluster_var cell_type --cluster $my_cell_type --group_var diagnosis --reference_group Control --alternative_group MS --covars Capbatch,age,sex --out_dir /home/jgagnon1/scRNAseq/downsample --out_prefix MS_Nature_downsample_0R1" | qsub -N downsample0R1 -q cpu.q -l h_rt=48:00:00 -M jake.gagnon@biogen.com -pe thread 8 -cwd -j yes -m ea

done


	

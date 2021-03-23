#!/bin/bash

cell_names="Glia-MIX IN-SST EN-L4 OL Microglia IN-PV IN-SV2C EN-PYR Astrocytes OPC IN-VIP EN-L2-3 EN-MIX EN-L5-6"
#cell_names="OL"
#cell_names="Glia-MIX"

for my_cell_type in $cell_names
do
    for my_split in $(seq 16)
    do
	echo "source /etc/profile.d/modules_bash.sh; module load R/3.5.1; Rscript --vanilla glmm300k_eugenia.R --cluster $my_cell_type --in_dir /home/jgagnon1/jupyter-notebook-dir/CALvsNWM_meta_filtered/ --out_dir /home/jgagnon1/scRNAseq/300k_results_CALvsNWM --split $my_split" |qsub -N run_300k_CALNWM -q cpu.q -l h_rt=48:00:00 -M jake.gagnon@biogen.com -cwd -j yes -m ea
    done
done

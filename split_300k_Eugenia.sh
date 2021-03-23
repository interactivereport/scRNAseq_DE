#!/bin/bash
#$ -N split_300k_Eugenia
#$ -q cpu.q
#$ -l h_rt=4:00:00
#$ -cwd
#$ -j yes
#$ -M jake.gagnon@biogen.com
#$ -m ea

source /etc/profile.d/modules_bash.sh 
module load R/3.5.1
Rscript --vanilla split_300k_eugenia.R 
exit 0

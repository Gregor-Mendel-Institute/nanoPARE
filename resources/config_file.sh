#!/bin/bash

root_dir=/lustre/scratch/users/michael.schon/EndGraph
bash_dir=/lustre/scratch/users/michael.schon/EndGraph/scripts/bash_scripts
python_dir=/lustre/scratch/users/michael.schon/EndGraph/scripts/python_scripts
r_dir=/lustre/scratch/users/michael.schon/EndGraph/scripts/r_scripts
resource_dir=/lustre/scratch/users/michael.schon/EndGraph/resources
temp_dir=/lustre/scratch/users/michael.schon/EndGraph/temp
log_dir=/lustre/scratch/users/michael.schon/EndGraph/log
results_dir=/lustre/scratch/users/michael.schon/EndGraph/results

genome_fasta=/lustre/scratch/users/michael.schon/EndGraph/resources/genome.fasta
annotation_gff=/lustre/scratch/users/michael.schon/EndGraph/resources/annotation.gff
TSS_PLUS=/lustre/scratch/users/michael.schon/EndGraph/resources/TSS_plus.bedgraph
TSS_MINUS=/lustre/scratch/users/michael.schon/EndGraph/resources/TSS_minus.bedgraph
TES_PLUS=/lustre/scratch/users/michael.schon/EndGraph/resources/TES_plus.bedgraph
TES_MINUS=/lustre/scratch/users/michael.schon/EndGraph/resources/TES_minus.bedgraph
BODY_PLUS=/lustre/scratch/users/michael.schon/EndGraph/resources/BODY_plus.bedgraph
BODY_MINUS=/lustre/scratch/users/michael.schon/EndGraph/resources/BODY_minus.bedgraph

LMOD=0
SETUP=true
ITERATIONS=5
KERNEL=laplace



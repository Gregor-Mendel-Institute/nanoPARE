#!/bin/bash

root_dir=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph
bash_dir=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/scripts/bash_scripts
python_dir=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/scripts/python_scripts
r_dir=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/scripts/r_scripts
resource_dir=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/resources
temp_dir=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/temp
log_dir=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/log
results_dir=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/results

genome_fasta=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/resources/genome.fasta
annotation_gff=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/resources/annotation.gff
TSS_PLUS=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/resources/TSS_plus.bedgraph
TSS_MINUS=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/resources/TSS_minus.bedgraph
TES_PLUS=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/resources/TES_plus.bedgraph
TES_MINUS=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/resources/TES_minus.bedgraph
BODY_PLUS=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/resources/BODY_plus.bedgraph
BODY_MINUS=/home/GMI/michael.schon/Desktop/lab/members/Schon/Writing/nanoPARE/EndGraph/resources/BODY_minus.bedgraph

LMOD=0
SETUP=false
ITERATIONS=5
KERNEL=laplace



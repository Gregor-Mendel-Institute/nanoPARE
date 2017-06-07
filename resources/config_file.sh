#!/bin/bash

root_dir=/lustre/scratch/users/michael.schon/bookend
bash_dir=/lustre/scratch/users/michael.schon/bookend/scripts/bash_scripts
python_dir=/lustre/scratch/users/michael.schon/bookend/scripts/python_scripts
r_dir=/lustre/scratch/users/michael.schon/bookend/scripts/r_scripts
resource_dir=/lustre/scratch/users/michael.schon/bookend/resources
temp_dir=/lustre/scratch/users/michael.schon/bookend/temp
log_dir=/lustre/scratch/users/michael.schon/bookend/log
results_dir=/lustre/scratch/users/michael.schon/bookend/results

fastq_dir=/lustre/scratch/users/michael.schon/bookend/resources/fastq_files
reference_table=/lustre/scratch/users/michael.schon/bookend/resources/reference.table
genome_fasta=/lustre/scratch/users/michael.schon/bookend/resources/genome.fasta
annotation_gff=/lustre/scratch/users/michael.schon/bookend/resources/annotation.gff

LMOD=1
CPUS=1
RAM=30
ISPCR=AAGCAGTGGTATCAACGCAGAGTAC
TN5_1=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
TN5_2=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

PBS_P=rnaseq_nod
PBS_q=workq



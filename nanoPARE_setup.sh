#!/bin/bash

function usage() {
cat <<HELP

### nanoPARE Setup ###
Usage: ./nanoPARE_setup.sh [options]

Sets up a collection of files required to run the nanoPARE tools:
 1) Writes a length.table recording the length of each chromosome in the genome FASTA
 2) Generates a STAR genome index using the reference genome and transcriptome
 3) Writes a BED file of putative TSO strand invasion sites based on mask_sequences.table  
 4) Writes a transcriptome FASTA file by combinding -G and -A  
 5) Parses reference annotations to identify 5'-terminal exons and single-exon transcripts  

The genome index is written to temp/genome_index.
Other setup files are added to the resources/ directory.

Optional arguments:
-R | --reference     Reference table (default: resources/reference.table)
-G | --genome        Genome FASTA file (default: resources/genome.fasta)
-A | --annotation    Transcript GFF file (default: resources/annotation.gff)
--lmod               Load required modules with Lmod (default: false)
--ram                Amount of available RAM in gigabytes (default: 30)
--cpus               Number of cores available for multithreaded programs (default: 1)

All steps of this pipeline access the default files for -R, -G, and -A, respectively.
You can replace these 3 items in the resources folder for simplicity.

HELP
}


################
# CONFIG SETUP #
################
# Storing all default global environment variables
if [ -z "$root_dir" ]
then
    root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" # Temporary $root_dir relative to this script's location if it isn't already in the environment
fi
bash_dir=$root_dir/scripts/bash_scripts
python_dir=$root_dir/scripts/python_scripts
resource_dir=$root_dir/resources
temp_dir=$root_dir/temp
log_dir=$root_dir/log
results_dir=$root_dir/results
genome_index_dir=$temp_dir/genome_index

. $bash_dir/read_cmdline.sh

# Set defaults for variables if they are not already in the environment
if [ -z "$GENOME_FASTA" ]
then
    GENOME_FASTA=$resource_dir/genome.fasta # If not already in environment, set as default value
fi
if [ -z "$REFERENCE_TABLE" ]
then
    REFERENCE_TABLE=$resource_dir/reference.table
fi
if [ -z "$ANNOTATION_GFF" ]
then
    ANNOTATION_GFF=$resource_dir/annotation.gff
fi
if [ -z "$LMOD" ]
then
    LMOD=0
fi
if [ -z "$CPUS" ]
then
    CPUS=1
fi

. $bash_dir/list_settings.sh

# Environment modules to load with Lmod (if option --lmod is passed)
REQUIRED_MODULES=( --rna-star --bedtools --python2018 )
. $bash_dir/load_modules.sh

cd $root_dir
echo "####################################"
echo "### SETUP: GENERATE GENOME INDEX ###"
echo "####################################"
# A table of chromosome names and lengths derived from the input FASTA genome file
length_table_command="python $python_dir/fasta_lengths.py $GENOME_FASTA > $resource_dir/length.table"
echo "$length_table_command"
eval "$length_table_command"
echo "Length table generated."
echo " "

length_array=($(cut -f 2 $resource_dir/length.table | tr '\n' ' '))
genome_length=$(echo ${length_array[@]} | tr ' ' '+' | bc)
index_base_number=$(echo "l($genome_length)/l(2)/2 -1" | bc -l | cut -d '.' -f 1)
if [ $index_base_number -gt 14 ]
then
    index_base_number=14
fi

# Genome index required by STAR for short read alignments
rm -rf $genome_index_dir
mkdir -p $genome_index_dir

genome_index_command="STAR \
--runThreadN $CPUS \
--runMode genomeGenerate \
--genomeDir $genome_index_dir \
--outFileNamePrefix $log_dir/ \
--genomeFastaFiles $GENOME_FASTA \
--genomeSAindexNbases $index_base_number \
--sjdbGTFfile $ANNOTATION_GFF"

echo "$genome_index_command"
eval "$genome_index_command"
echo "Genome index complete."
echo " " 

echo "### GENERATE TSO MASKING FILES ###"

python $python_dir/fasta_sequence_search.py \
    $GENOME_FASTA \
    $resource_dir/mask_sequences.table \
    -O $resource_dir
echo "Masking BED files generated."
echo " "

echo "Generating transcriptome FASTA file."
python $python_dir/gtf_to_fasta.py \
    -G $GENOME_FASTA \
    -A $ANNOTATION_GFF \
    > $resource_dir/transcriptome.fasta

echo "### GENERATE ANNOTATION CLASS REFERENCE FILES ###"

cd $resource_dir
echo "   Getting transcript-level exons"
grep -P "exon\t" $ANNOTATION_GFF |\
    sed 's/\t[^\t]*transcript_id=\([^;]*\);.*/\t\1\t/' |\
    sed 's/Parent=\(transcript:\)\?\([^;]*\).*/\2/' |\
    sort -k 9 > class.exons_by_transcript.gff

echo "   Getting 5'-most exons"
python $python_dir/bed_deduplicate.py\
    -F 8 -S upstream --startline 3 --endline 4 --strandline 6 class.exons_by_transcript.gff > class.terminal_exons_by_transcript.gff

sed 's/\.[0-9]\{1,\}$//' class.terminal_exons_by_transcript.gff | sort -k 9 > class.terminal_exons_by_gene.gff
echo "   Converting transcript-level to gene-level annotations"
grep -P "exon\t" $ANNOTATION_GFF | sed 's/gene_id=\([^;]*\);.*/\1\t/' | sed 's/Parent=\(transcript:\)\?\([^;\.]*\).*/\2/' | sort -k 9 > class.exons_by_gene.gff

echo "   Reformatting to GFF files to BED files"
bedtools sort -i class.exons_by_gene.gff |\
    awk '{ printf $1"\t"$4-1"\t"$5"\t"$9"\t0\t"$7"\n" }' |\
    bedtools merge -s -c 4 -o distinct -delim ";" |\
    awk '{ printf $1"\t"$2"\t"$3"\t"$5"\t0\t"$4"\n" }' |\
    sed 's/\(.\+\)\t\([^;]*\?\);.\+\t/\1\t\2\t/' |\
    bedtools sort | sort -k 4 > class.exons_by_gene.geneorder.bed

bedtools sort -i class.exons_by_gene.geneorder.bed > class.exons_by_gene.bed

echo "   Reformatting to GFF files to BED files"
bedtools sort -i class.terminal_exons_by_gene.gff |\
    awk '{ printf $1"\t"$4-1"\t"$5"\t"$9"\t0\t"$7"\n" }' |\
    bedtools merge -s -c 4 -o distinct -delim ";" |\
    awk '{ printf $1"\t"$2"\t"$3"\t"$5"\t0\t"$4"\n" }' |\
    sed 's/\(.\+\)\t\([^;]*\?\);.\+\t/\1\t\2\t/' |\
    bedtools sort > class.terminal_exons_by_gene.bed

python $python_dir/bed_deduplicate.py \
    --only_unique class.exons_by_gene.geneorder.bed -F 3 |\
    bedtools sort > class.single_exon_genes.bed

echo "   Cleaning up temporary files"
rm class.terminal_exons_by_gene.gff class.exons_by_gene.gff class.exons_by_gene.geneorder.bed
rm class.terminal_exons_by_transcript.gff class.exons_by_transcript.gff

echo "Setup complete."
exit 0




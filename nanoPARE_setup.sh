#!/bin/bash

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
r_dir=$root_dir/scripts/r_scripts
resource_dir=$root_dir/resources
temp_dir=$root_dir/temp
log_dir=$root_dir/log
results_dir=$root_dir/results
genome_index_dir=$temp_dir/genome_index
ENDMAP_DIR=$results_dir/EndMap
ENDGRAPH_DIR=$results_dir/EndGraph

# Set defaults for variables if they are not already in the environment
if [ -z "$JOB_NUMBER" ]
then
    JOB_NUMBER=${PBS_ARRAY_INDEX} # Imports the PBS job array number if it exists. Can be overridden with the commandline argument -J $JOB_NUMBER
fi
if [ -z "$genome_fasta" ]
then
    genome_fasta=$resource_dir/genome.fasta # If not passed already in environment, set as default value
fi
if [ -z "$reference_table" ]
then
    reference_table=$resource_dir/reference_table_EndGraph.txt
fi
if [ -z "$endmap_reference_table" ]
then
    endmap_reference_table=$resource_dir/reference_table_EndMap.txt
fi

if [ -z "$annotation_gff" ]
then
    annotation_gff=$resource_dir/annotation.gff
fi
if [ -z "$annotation_subset" ]
then
    annotation_subset=$resource_dir/annotation_subset.txt
fi
if [ -z "$sample_name" ]
then
    sample_name="sample"
fi
if [ -z "$LMOD" ]
then
    LMOD=0
fi
if [ -z "$CPUS" ]
then
    CPUS=1
fi
if [ -z "$SETUP" ]
then
    SETUP=false
fi

KERNEL='laplace'
WEIGHTED='true'

############################
# READING THE COMMAND LINE #
############################
# Taking the default variables above, modifying them with the commandline 
# arguments (see read_cmdline.sh), and writing a config file

. $bash_dir/read_cmdline.sh

echo "################"
echo "### ENDGRAPH ###"
echo "################"
echo " "
echo "Config settings:"
. $bash_dir/list_settings.sh

# Environment modules to load with Lmod (if option --lmod is passed)
REQUIRED_MODULES=( --bedtools --python )
. $bash_dir/load_modules.sh
echo " "

input_mapper=$(sed -n "$JOB_NUMBER"p $reference_table) #read mapping file
input_array=($input_mapper)

line_number=${input_array[0]}  # Line number of reference table (should match job number)
sample_name=${input_array[1]}  # Name of the positive sample
body_name=${input_array[2]}    # Name of the background sample
library_type=${input_array[3]} # Options: 5P, 3P. Type of the positive sample

sample_dir=$ENDGRAPH_DIR/$sample_name
mkdir -p $sample_dir


echo "#######################################"
echo "### SETUP: GENERATE ANNOTATION INFO ###"
echo "#######################################"
echo " "
# A table of chromosome names and lengths derived from the input FASTA genome file
length_table_command="python $python_dir/fasta_lengths.py $genome_fasta > $resource_dir/length.table"
echo "$length_table_command"
eval "$length_table_command"
echo "Length table generated."

python $python_dir/fasta_sequence_search.py \
    $genome_fasta \
    $resource_dir/mask_sequences.table \
    -O $resource_dir
echo "Masking BED files generated."

echo "Setting up exon reference..."
cd $resource_dir
echo "   Getting transcript-level exons"
grep -P "exon\t" $annotation_gff |\
    sed 's/\t[^\t]*transcript_id=\([^;]*\);.*/\t\1\t/' |\
    sed 's/Parent=\(transcript:\)\?\([^;]*\).*/\2/' |\
    sort -k 9 > exons_by_transcript.gff

echo "   Getting transcript-level CDS"
grep -P "CDS\t" $annotation_gff |\
    sed 's/\t[^\t]*transcript_id=\([^;]*\);.*/\t\1\t/' |\
    sed 's/\(Parent\|ID\)=\(CDS:\)\?\([^;\.]*\).*/\3/' |\
    sort -k 9 > CDS_by_transcript.gff
echo "   Getting 5'-most exons"
python $python_dir/bed_deduplicate.py\
    -F 8 -S upstream --startline 3 --endline 4 --strandline 6 exons_by_transcript.gff > terminal_exons_by_transcript.gff

sed 's/\.[0-9]\{1,\}$//' terminal_exons_by_transcript.gff | sort -k 9 > terminal_exons_by_gene.gff
echo "   Converting transcript-level to gene-level annotations"
grep -P "exon\t" $annotation_gff | sed 's/gene_id=\([^;]*\);.*/\1\t/' | sed 's/Parent=\(transcript:\)\?\([^;\.]*\).*/\2/' | sort -k 9 > exons_by_gene.gff
grep -P "CDS\t" $annotation_gff | sed 's/gene_id=\([^;]*\);.*/\1\t/' | sed 's/\(Parent\|ID\)=\(CDS:\)\?\([^;\.]*\).*/\3/' | sort -k 9 > CDS_by_gene.gff

echo "   Reformatting to BED file"
bedtools sort -i exons_by_gene.gff |\
        awk '{ printf $1"\t"$4-1"\t"$5"\t"$9"\t"$9"\t"$7"\n" }' |\
        bedtools merge -s -nms |\
        sed 's/\(.\+\)\t\([^;]*\?\);.\+\t/\1\t\2\t/' |\
        awk '{ printf $1"\t"$2"\t"$3"\t"$4"\t0\t"$5"\n" }' |\
        bedtools sort | sort -k 4 > exons_by_gene.bed

bedtools sort -i CDS_by_gene.gff |\
        awk '{ printf $1"\t"$4-1"\t"$5"\t"$9"\t"$9"\t"$7"\n" }' |\
        bedtools merge -s -nms |\
        sed 's/\(.\+\)\t\([^;]*\?\);.\+\t/\1\t\2\t/' |\
        awk '{ printf $1"\t"$2"\t"$3"\t"$4"\t0\t"$5"\n" }' |\
        bedtools sort | sort -k 4 > CDS_by_gene.bed

echo "   Getting genes with single-exon transcripts"
python $python_dir/bed_deduplicate.py \
    --only_unique exons_by_gene.bed -F 3 > single_exon_genes.bed


echo "Setup complete."
exit 0




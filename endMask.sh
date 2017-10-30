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
endmap_dir=$results_dir/EndMap
endgraph_dir=$results_dir/EndGraph
genome_index_dir=$temp_dir/genome_index

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
    reference_table=$resource_dir/reference_table_EndMask.txt
fi
if [ -z "$annotation_gff" ]
then
    annotation_gff=$resource_dir/annotation.gff
fi
if [ -z "$annotation_subset" ]
then
    annotation_subset=$resource_dir/annotation_subset.txt
fi
if [ -z "$SAMPLE_NAME" ]
then
    SAMPLE_NAME="sample"
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
if [ -z "$output_folder" ]
then
    output_folder=$results_dir/EndMask
fi

############################
# READING THE COMMAND LINE #
############################
# Taking the default variables above, modifying them with the commandline 
# arguments (see read_cmdline.sh), and writing a config file

. $bash_dir/read_cmdline.sh

echo "###############"
echo "### ENDMASK ###"
echo "###############"
echo " "
echo "Config settings:"
. $bash_dir/list_settings.sh


# Environment modules to load with Lmod (if option --lmod is passed)
REQUIRED_MODULES=( --bt )
. $bash_dir/load_modules.sh
echo " "

input_mapper=$(sed -n "$JOB_NUMBER"p $reference_table) #read mapping file
input_array=($input_mapper)

line_number=${input_array[0]}  # Line number of reference table (should match job number)
SAMPLE_NAME=${input_array[1]}  # name of sample type
sample_list=$(echo ${input_array[2]} | tr ',' ' ')  # parse the comma-separated list of sample names

#####################
# ENVIRONMENT SETUP #
#####################

length_table=$resource_dir/length.table
mask=$python_dir/bedgraph_mask.py
coverage=$python_dir/bed_feature_coverage.py
samples=($sample_list)
endmap_folder=$results_dir/EndMap
endgraph_folder=$temp_dir
output_folder=$results_dir/EndMask
data_folder=$output_folder/$SAMPLE_NAME

echo "Samples: ${samples[@]}"
echo "Sample type: $SAMPLE_NAME"
mkdir -p $data_folder
cd $data_folder

echo "Setting up exon reference..."
grep -P "\texon\t" $annotation_gff > reference_exons.gff
sed 's/gene_id=\([^;]*\);.*/\1\t/' reference_exons.gff > exons_by_gene.gff
awk '{printf $1"\t"$4-1"\t"$5"\t"$9"\t"$6"\t"$7"\t"$8"\n"}' exons_by_gene.gff > exons.bed

##################
# MERGE FEATURES #
##################

echo "Merging capped and uncapped feature files..."
for captype in capped uncapped
do
    rm -f $SAMPLE_NAME.$captype.bed
    touch $SAMPLE_NAME.$captype.bed
    for s in ${samples[@]}
    do
        sed "s/\tTSS\./\t"$s"\.TSS\./" \
            $endgraph_folder/$s/"$captype"_features.bed \
            >> $SAMPLE_NAME.$captype.bed
    done

    bedtools merge \
        -s \
        -nms \
        -n \
        -d 0 \
        -i $SAMPLE_NAME.$captype.bed \
        > "$SAMPLE_NAME"_merged.$captype.bed

    grep -v -P '\t1\t' "$SAMPLE_NAME"_merged.$captype.bed \
        > $SAMPLE_NAME.rep.$captype.bed

    rm $SAMPLE_NAME.$captype.bed "$SAMPLE_NAME"_merged.$captype.bed

    bedfile=$SAMPLE_NAME.rep.$captype.bed

    bedtools closest \
        -s \
        -D b \
        -t first \
        -a $bedfile \
        -b reference_exons.gff \
        > "$SAMPLE_NAME"_closest_"$captype".bed

    sed 's/gene_id=\([^;]*\);.*\t/\1\t/' "$SAMPLE_NAME"_closest_"$captype".bed > "$SAMPLE_NAME"_"$captype"_closest_gene.bed

    awk '{printf $1"\t"$2"\t"$3"\t"$15"\t"$16"\t"$6"\t"$5"\n"}' "$SAMPLE_NAME"_"$captype"_closest_gene.bed \
        > "$SAMPLE_NAME"_gene_"$captype".bed
    
    awk -F'[\t]' \
        'function abs(v) {return v < 0 ? -v : v};
        {if (abs($5) <= 50){ print }}' \
        "$SAMPLE_NAME"_gene_"$captype".bed \
        > "$SAMPLE_NAME"_"$captype"_overlap.bed

    rm -f "$SAMPLE_NAME"_closest_"$captype".bed "$SAMPLE_NAME"_"$captype"_closest_gene.bed
done

echo "Cap masking bedgraph files..."

capped_bedgraphs=()
uncapped_bedgraphs=()
names=()

for s in ${samples[@]}
do
    python $mask \
    -P $endgraph_folder/$s/TSS_plus_mask.bedgraph \
    -M $endgraph_folder/$s/TSS_minus_mask.bedgraph \
    -PO $data_folder/"$s"_capmask_plus.bedgraph \
    -MO $data_folder/"$s"_capmask_minus.bedgraph \
    -I $SAMPLE_NAME.rep.capped.bed \
    -L $length_table
    
    names+=( $s.plus $s.minus )
    capped_bedgraphs+=( $endgraph_folder/$s/TSS_plus_mask.bedgraph $endgraph_folder/$s/TSS_minus_mask.bedgraph )
    uncapped_bedgraphs+=( $data_folder/"$s"_capmask_plus.bedgraph $data_folder/"$s"_capmask_minus.bedgraph )
    
done

###########################
# CALCULATE READ COVERAGE #
###########################

echo "Making gene-level bed file..."
cat "$SAMPLE_NAME"_capped_overlap.bed "$SAMPLE_NAME"_uncapped_overlap.bed exons.bed > "$SAMPLE_NAME"_gene_features.bed

echo "Calculating gene-level total and uncapped read coverage..."
python $coverage \
    -F "$SAMPLE_NAME"_capped_overlap.bed \
    -I ${capped_bedgraphs[@]} \
    -N ${names[@]} \
    -O "$SAMPLE_NAME"_capped_coverage.tsv \
    -L $length_table 

python $coverage \
    -F "$SAMPLE_NAME"_gene_features.bed \
    -I ${capped_bedgraphs[@]} \
    -N ${names[@]} \
    -O "$SAMPLE_NAME"_total_coverage.tsv \
    -L $length_table \
    -G $genome_fasta \
    --g_content "$SAMPLE_NAME"_g_content.tsv

python $coverage \
    -F "$SAMPLE_NAME"_gene_features.bed \
    -I ${uncapped_bedgraphs[@]} \
    -N ${names[@]} \
    -O "$SAMPLE_NAME"_uncapped_coverage.tsv \
    -L $length_table




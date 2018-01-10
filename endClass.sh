#!/bin/bash

# EndClass
# Classification of 5-prime end features as capped or noncapped
# based on the presence or absence of upstream untemplated guanosines 
# 1) Filters 5P features for an experiment for those replicable in >=2 expts
# 2) Counts the proportion of reads in each 5P feature with upstream untemplated G (uuG)
# 3) Splits 5P features into "capped" (>=10% uuG) and "noncapped" (<10% uuG)

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
    reference_table=$resource_dir/reference_table_EndClass.txt
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
if [ -z "$output_folder" ]
then
    output_folder=$results_dir/EndClass
fi

############################
# READING THE COMMAND LINE #
############################
# Taking the default variables above, modifying them with the commandline 
# arguments (see read_cmdline.sh), and writing a config file

. $bash_dir/read_cmdline.sh

echo "###############"
echo "### ENDCLASS ###"
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
sample_name=${input_array[1]}  # name of sample type
sample_list=$(echo ${input_array[2]} | tr ',' ' ')  # parse the comma-separated list of sample names

#####################
# ENVIRONMENT SETUP #
#####################

length_table=$resource_dir/length.table
mask=$python_dir/bedgraph_mask.py
coverage=$python_dir/bed_feature_coverage.py
samples=($sample_list)
endmap_folder=$results_dir/EndMap
data_folder=$output_folder/$sample_name

echo "Samples: ${samples[@]}"
echo "Sample type: $sample_name"

mkdir -p $data_folder
cd $data_folder

echo "Setting up exon reference..."
grep -P "exon\t" $annotation_gff | sed 's/gene_id=\([^;]*\);.*/\1\t/' | sed 's/Parent=\([^;\.]*\).*/\1/' > exons_by_gene.gff
awk '{printf $1"\t"$4-1"\t"$5"\t"$9"\t"$6"\t"$7"\t"$8"\n"}' exons_by_gene.gff > exons.bed

##################
# MERGE FEATURES #
##################

echo "Merging feature files..."
rm -f $sample_name.all.bed
touch $sample_name.all.bed
plus_list=()
minus_list=()
plus_list_uug=()
minus_list_uug=()
for s in ${samples[@]}
do
    sed "s/\t5P\./\t"$s"\.5P\./" \
        $temp_dir/$s/$s.end_features.bed \
        >> $sample_name.all.bed
    
    plus_list+=($temp_dir/$s/5P_plus_mask.bedgraph)
    minus_list+=($temp_dir/$s/5P_minus_mask.bedgraph)
    plus_list_uug+=($temp_dir/$s/uuG_plus_mask.bedgraph)
    minus_list_uug+=($temp_dir/$s/uuG_minus_mask.bedgraph)
done

echo "Bedgraph files +: "${plus_list[@]}
echo "Bedgraph files -: "${minus_list[@]}
python $python_dir/bedgraph_combine.py \
    -i ${plus_list[@]} \
    -o $sample_name.plus.bedgraph

python $python_dir/bedgraph_combine.py \
    -i ${minus_list[@]} \
    -o $sample_name.minus.bedgraph

echo "Bedgraph uuG files +: "${plus_list_uug[@]}
echo "Bedgraph uuG files -: "${minus_list_uug[@]}
python $python_dir/bedgraph_combine.py \
    -i ${plus_list_uug[@]} \
    -o $sample_name.uuG.plus.bedgraph

python $python_dir/bedgraph_combine.py \
    -i ${minus_list_uug[@]} \
    -o $sample_name.uuG.minus.bedgraph

echo "Merged coverage files generated."

# Merge all touching/overlapping features from all reps
bedtools merge \
    -s \
    -nms \
    -n \
    -d 0 \
    -i $sample_name.all.bed \
    > $sample_name.merged.bed

# Keep only replicable features (present in more than one library)
grep -v -P '\t1\t' $sample_name.merged.bed \
    > $sample_name.rep.bed


rm $sample_name.all.bed $sample_name.merged.bed

bedfile=$sample_name.rep.bed

# Locate the nearest (sense) gene in exons_by_gene.gff
bedtools closest \
    -s \
    -D b \
    -t first \
    -a $bedfile \
    -b exons_by_gene.gff \
    > $sample_name.closest_gene.bed

# Ath_chr1        310073  310434  Ath_chr1.636;Ath_chr1.672;Ath_chr1.652  3       +       Ath_chr1        TAIR10  exon    310316  310981  .       +       .       AT1G01900       0
awk '{printf $1"\t"$2"\t"$3"\t5P."NR"\t"$16"\t"$6"\t"$15"\n"}' $sample_name.closest_gene.bed \
    > $sample_name.gene.bed

# Adding a merged peak position ot each replicable feature, keeping both the closest gene ID and distance
python $python_dir/bed_find_peaks.py \
    -I $sample_name.gene.bed \
    -O $sample_name.all_features.bed \
    -BP $sample_name.plus.bedgraph \
    -BM $sample_name.minus.bedgraph \
    -L $length_table \
    -V pass

rm -f $sample_name.closest_gene.bed $sample_name.gene.bed

echo "Splitting capped and noncapped features..."

python $python_dir/bed_uug_filter.py \
    -C $sample_name.capped.bed \
    -U $sample_name.noncapped.bed \
    --min_uug 0.1 \
    $sample_name.all_features.bed \
    $sample_name.plus.bedgraph \
    $sample_name.minus.bedgraph \
    $sample_name.uuG.plus.bedgraph \
    $sample_name.uuG.minus.bedgraph

echo "Cap masking bedgraph files..."

capped_bedgraphs=()
noncapped_bedgraphs=()
names=()

for s in ${samples[@]}
do
    python $mask \
    -P $temp_dir/$s/5P_plus_mask.bedgraph \
    -M $temp_dir/$s/5P_minus_mask.bedgraph \
    -PO $data_folder/$s.capmask.plus.bedgraph \
    -MO $data_folder/$s.capmask.minus.bedgraph \
    -I $sample_name.capped.bed \
    -L $length_table
    
    names+=( $s.plus $s.minus )
    capped_bedgraphs+=( $temp_dir/$s/5P_plus_mask.bedgraph $temp_dir/$s/5P_minus_mask.bedgraph )
    noncapped_bedgraphs+=( $data_folder/$s.capmask.plus.bedgraph $data_folder/$s.capmask.minus.bedgraph )
done

echo "Cap masking merged bedgraph..."
python $mask \
-P $sample_name.plus.bedgraph \
-M $sample_name.minus.bedgraph \
-PO $sample_name.capmask.plus.bedgraph \
-MO $sample_name.capmask.minus.bedgraph \
-I $sample_name.capped.bed \
-L $length_table

###########################
# CALCULATE READ COVERAGE #
###########################

echo "Making gene-level bed file..."
# Subset 5' end features for only those that are
# (1) capped, and
# (2) within 50nt of an existing gene annotation

awk -F'[\t]' \
    'function abs(v) {return v < 0 ? -v : v};
    {if (abs($5) <= 50){ print }}' \
    $sample_name.capped.bed \
    > $sample_name.overlapping.capped.bed

cat $sample_name.overlapping.capped.bed exons.bed > "$sample_name"_gene_features.bed
rm $sample_name.overlapping.capped.bed

echo "Calculating gene-level total and uncapped read coverage..."
python $coverage \
    -F $sample_name.capped.bed \
    -I ${capped_bedgraphs[@]} \
    -N ${names[@]} \
    -O "$sample_name"_capped_coverage.tsv \
    -L $length_table 

python $coverage \
    -F $sample_name.noncapped.bed \
    -I ${capped_bedgraphs[@]} \
    -N ${names[@]} \
    -O "$sample_name"_noncapped_coverage.tsv \
    -L $length_table

python $coverage \
    -F "$sample_name"_gene_features.bed \
    -I ${capped_bedgraphs[@]} \
    -N ${names[@]} \
    -O "$sample_name"_total_coverage.tsv \
    -L $length_table

python $coverage \
    -F "$sample_name"_gene_features.bed \
    -I ${noncapped_bedgraphs[@]} \
    -N ${names[@]} \
    -O "$sample_name"_capmasked_coverage.tsv \
    -L $length_table



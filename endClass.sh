#!/bin/bash

function usage() {
cat <<HELP

### EndClass ###
Usage: ./endClass.sh [options] -T|--type <sample_type>

Labeling and classification of 5'-end features as capped or noncapped
 1) Labels all 5P features based on their relationship to the nearest transcriptome annotation
 2) Filters 5P features for only those replicable in >=2 experiments
 3) Counts the proportion of reads in each 5P feature with upstream untemplated G (uuG)
 4) Splits 5P features into "capped" (>=10% uuG) and "noncapped" (<10% uuG)

Expects EndMap and EndGraph to be run before and their results to be in directories
/results/EndMap/sample_name
and
/results/EndGraph/sample_name
respectively, for all sample_names of the chosen sample_type that appear in reference.table

Optional arguments:
-R | --reference     Reference table (default: resources/reference.table)
-G | --genome        Genome FASTA file (default: resources/genome.fasta)
-A | --annotation    Transcript GFF file (default: resources/annotation.gff)
--lmod               Load required modules with Lmod (default: false)
--cpus               Number of cores available for multithreaded programs (default: 1)
--uug                Minimum proportion uuG required to classify a feature as capped (default: 0.1)
--upstream           Maximum distance upstream (in nucleotides) to associate a 5P feature (default: 500)

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
endmap_dir=$results_dir/EndMap
endgraph_dir=$results_dir/EndGraph
endclass_dir=$results_dir/EndClass
genome_index_dir=$temp_dir/genome_index

############################
# READING THE COMMAND LINE #
############################
# Taking the default environment above, add the commandline
# arguments (see read_cmdline.sh), and write a config file
if [ $# -eq 0 ]; then
    if [ -z $SAMPLE_ARRAY ]
    then
        usage
        exit 1
    fi
else
    . $bash_dir/read_cmdline.sh
fi

# Set defaults for variables if they are not already in the environment
if [ -z "$JOB_NUMBER" ]
then
    JOB_NUMBER=$(echo "${PBS_ARRAY_INDEX}-1" | bc) # Imports the PBS job array number if it exists. Can be overridden with the commandline argument -J $JOB_NUMBER
fi
if [ -z "$SAMPLE_TYPE" ]
then
    if [ -z "$SAMPLE_ARRAY" ]
    then
        echo "ERROR: Please input a sample type (-T|--type)"
        exit 1
    fi
    SAMPLES=( $(echo $SAMPLE_ARRAY | tr '!' ' ' ) )
    SAMPLE_TYPE=${SAMPLES[$JOB_NUMBER]}
fi
if [ -z "$GENOME_FASTA" ]
then
    GENOME_FASTA=$resource_dir/genome.fasta # If not passed already in environment, set as default value
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
if [ -z "$UUG" ]
then
    UUG=0.1
fi
if [ -z "$UPSTREAM" ]
then
    UPSTREAM=500
fi

echo "Config settings:"
. $bash_dir/list_settings.sh

# Environment modules to load with Lmod (if option --lmod is passed)
REQUIRED_MODULES=( --bedtools )
. $bash_dir/load_modules.sh
echo " "


#####################
# ENVIRONMENT SETUP #
#####################

sample_name_array=( $(cut -d ' ' -f 4 $REFERENCE_TABLE) )
sample_type_array=( $(cut -d ' ' -f 5 $REFERENCE_TABLE) )
picked_from_table=()
for index in ${!sample_type_array[*]}
do
    type=${sample_type_array[$index]}
    if [ $type == $SAMPLE_TYPE ]
    then
        name=${sample_name_array[$index]}
        picked_from_table+=( $name )
    fi
done
samples=( $(echo "${picked_from_table[@]}" | tr [:blank:] \\n | sort -u) )

length_table=$resource_dir/length.table
mask=$python_dir/bedgraph_mask.py
coverage=$python_dir/bed_feature_coverage.py
data_folder=$endclass_dir/$SAMPLE_TYPE

echo "Samples: ${samples[@]}"
echo "Sample type: $SAMPLE_TYPE"

mkdir -p $data_folder
cd $data_folder

##################
# MERGE FEATURES #
##################
echo "Merging feature files..."
rm -f $SAMPLE_TYPE.all.bed
touch $SAMPLE_TYPE.all.bed
plus_list=()
plus_names=()
minus_list=()
minus_names=()
plus_list_uug=()
minus_list_uug=()
for s in ${samples[@]}
do
    sed "s/\t5P\./\t"$s"\.5P\./" \
        $endgraph_dir/$s/$s.end_features.bed \
        >> $SAMPLE_TYPE.all.bed
    
    plus_list+=($endgraph_dir/$s/5P.plus_mask.bedgraph)
    plus_names+=($s.plus)
    minus_list+=($endgraph_dir/$s/5P.minus_mask.bedgraph)
    minus_names+=($s.minus)
    plus_list_uug+=($endgraph_dir/$s/uG.plus_mask.bedgraph)
    minus_list_uug+=($endgraph_dir/$s/uG.minus_mask.bedgraph)
done

echo "Bedgraph files +: "${plus_list[@]}
echo "Bedgraph files -: "${minus_list[@]}
python $python_dir/bedgraph_combine.py \
    -i ${plus_list[@]} \
    -o $SAMPLE_TYPE.plus.bedgraph

python $python_dir/bedgraph_combine.py \
    -i ${minus_list[@]} \
    -o $SAMPLE_TYPE.minus.bedgraph

echo "Bedgraph uuG files +: "${plus_list_uug[@]}
echo "Bedgraph uuG files -: "${minus_list_uug[@]}
python $python_dir/bedgraph_combine.py \
    -i ${plus_list_uug[@]} \
    -o $SAMPLE_TYPE.uG.plus.bedgraph

python $python_dir/bedgraph_combine.py \
    -i ${minus_list_uug[@]} \
    -o $SAMPLE_TYPE.uG.minus.bedgraph

echo "Merged coverage files generated."

# Merge all touching/overlapping features from all reps
# Keep only replicable features (present in more than one library)
bedtools sort -i $SAMPLE_TYPE.all.bed |\
    bedtools merge -s -c 5 -o count -d 0 |\
    grep -v -P '\t1$' |\
    awk -v nm="$SAMPLE_TYPE" '{ printf $1"\t"$2"\t"$3"\t"nm"."NR"\t"$5"\t"$4"\n" }' > $SAMPLE_TYPE.rep.bed

rm $SAMPLE_TYPE.all.bed

bedfile=$SAMPLE_TYPE.rep.bed

# Adding a merged peak position to each replicable feature, keeping both the closest gene ID and distance
echo " "
echo "Finding feature peak positions."
python $python_dir/bed_find_peaks.py \
    -I $bedfile \
    -O $SAMPLE_TYPE.peaks.bed \
    -BP $SAMPLE_TYPE.plus.bedgraph \
    -BM $SAMPLE_TYPE.minus.bedgraph \
    -L $length_table --trim 1 \
    -V pass

# Classify features against the annotated list of genes
# 1: Generate a bed file of 1-nt features representing only the peaks of each 5P feature
awk '{ printf $1"\t"$2+$5"\t"$2+$5+1"\t"$4"\t"$5"\t"$6"\t"$2"\t"$3"\t"$5"\n" }'\
    $SAMPLE_TYPE.peaks.bed | bedtools sort > $SAMPLE_TYPE.only_peaks.bed
bedfile=$SAMPLE_TYPE.only_peaks.bed

# 2: Find the nearest overlapping/downstream exon element in the single-exon transcripts
echo "Locating nearest annotation to each feature:"
echo "  overlapping single-exon transcripts"
bedtools closest -s -id -D b -t first \
    -a $bedfile \
    -b $resource_dir/class.single_exon_genes.bed > $SAMPLE_TYPE.single.tmp.bed
grep -P '\t0$' $SAMPLE_TYPE.single.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$14"\t7\n" }' > $SAMPLE_TYPE.single_exon_overlap.bed
rm $SAMPLE_TYPE.single.tmp.bed

# 3: Find the nearest overlapping/downstream exon element of 5'-terminal exons
echo "  overlapping 5'-terminal exons"
bedtools closest -s -id -D b -t first \
    -a $bedfile \
    -b $resource_dir/class.terminal_exons_by_gene.bed > $SAMPLE_TYPE.terminal.tmp.bed
grep -P '\t0$' $SAMPLE_TYPE.terminal.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$14"\t6\n" }' > $SAMPLE_TYPE.terminal_exon_overlap.bed
grep -v -P '\t0$' $SAMPLE_TYPE.terminal.tmp.bed | awk -v up="$UPSTREAM" '{ if ( $16 >= -up ) printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$16"\t4\n" }' > $SAMPLE_TYPE.upstream.bed
rm $SAMPLE_TYPE.terminal.tmp.bed

# 3: Find the nearest overlapping/downstream exon element of all transcripts
echo "  Between exons"
bedtools closest -s -id -D b -t first \
    -a $bedfile \
    -b $resource_dir/class.exons_by_gene.bed > $SAMPLE_TYPE.nearest_downstream.tmp.bed
grep -P '\t0$' $SAMPLE_TYPE.nearest_downstream.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$16"\t5\n" }' > $SAMPLE_TYPE.exon_overlap.bed
bedtools intersect -v -wa -a $SAMPLE_TYPE.exon_overlap.bed -b $SAMPLE_TYPE.terminal_exon_overlap.bed > $SAMPLE_TYPE.internal_exon_overlap.bed
grep -v -P '\t0$' $SAMPLE_TYPE.nearest_downstream.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$16"\t0\n" }' > $SAMPLE_TYPE.is_upstream_of.bed
bedtools intersect -v -wa -a $SAMPLE_TYPE.is_upstream_of.bed -b $SAMPLE_TYPE.upstream.bed | bedtools sort > $SAMPLE_TYPE.nearest_downstream.bed
rm $SAMPLE_TYPE.nearest_downstream.tmp.bed $SAMPLE_TYPE.exon_overlap.bed

# 4: Find the nearest overlapping/upstream exon element of all transcripts
echo "  Upstream of annotations"
bedtools closest -s -iu -D b -t first \
    -a $bedfile \
    -b $resource_dir/class.exons_by_gene.bed > $SAMPLE_TYPE.nearest_upstream.tmp.bed
grep -v -P '\t0$' $SAMPLE_TYPE.nearest_upstream.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$16"\t0\n" }' > $SAMPLE_TYPE.is_downstream_of.bed
bedtools intersect -v -wa -a $SAMPLE_TYPE.is_downstream_of.bed -b $SAMPLE_TYPE.upstream.bed | bedtools sort > $SAMPLE_TYPE.nearest_upstream.bed
rm $SAMPLE_TYPE.nearest_upstream.tmp.bed $SAMPLE_TYPE.is_downstream_of.bed

# 5: Isolate features whose peaks are both upstream and downstream of exons of the same gene
echo "  Contained within introns"
bedtools closest -a $SAMPLE_TYPE.nearest_upstream.bed -b $SAMPLE_TYPE.nearest_downstream.bed \
    | awk '{ if ( $18 == 0 && $7 == $16 ) printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$17"\t3\n" }' \
    > $SAMPLE_TYPE.intronic.bed

cat $SAMPLE_TYPE.terminal_exon_overlap.bed $SAMPLE_TYPE.upstream.bed $SAMPLE_TYPE.internal_exon_overlap.bed $SAMPLE_TYPE.intronic.bed \
    | bedtools sort > $SAMPLE_TYPE.all_sense.bed
bedtools intersect -v -s -wa -a $bedfile -b $SAMPLE_TYPE.all_sense.bed | bedtools sort > $SAMPLE_TYPE.residuals.bed

# 6: Find the nearest overlapping/upstream exon element in antisense for all features that had no sense match above
echo "  Antisense to existing annotations"
bedtools closest -S -id -D b -t first \
    -a $SAMPLE_TYPE.residuals.bed \
    -b $resource_dir/class.exons_by_gene.bed > $SAMPLE_TYPE.antisense_a.tmp.bed
grep -P '\t0$' $SAMPLE_TYPE.antisense_a.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$16"\t2\n" }' | bedtools sort > $SAMPLE_TYPE.antisense.bed
grep -v -P '\t0$' $SAMPLE_TYPE.antisense_a.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$16"\t0\n" }' | bedtools sort > $SAMPLE_TYPE.anti_upstream_of.bed

# 7: Find the nearest overlapping/downstream exon element in antisense for all features that had no sense match above
bedtools closest -S -iu -D b -t first \
    -a $SAMPLE_TYPE.residuals.bed \
    -b $resource_dir/class.exons_by_gene.bed > $SAMPLE_TYPE.antisense_d.tmp.bed
grep -v -P '\t0$' $SAMPLE_TYPE.antisense_d.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$16"\t0\n" }' | bedtools sort > $SAMPLE_TYPE.anti_downstream_of.bed

# 5: Isolate features whose peaks are both upstream and downstream of antisense exons of the same gene
echo "  Antisense intronic"
bedtools closest -a $SAMPLE_TYPE.anti_upstream_of.bed -b $SAMPLE_TYPE.anti_downstream_of.bed \
    | awk '{ if ( $18 == 0 && $7 == $16 ) printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$17"\t1\n" }' \
    > $SAMPLE_TYPE.anti_intronic.bed

rm $SAMPLE_TYPE.antisense_a.tmp.bed $SAMPLE_TYPE.antisense_d.tmp.bed

# Collapse the decision tree into a single BED file, with each feature
# represented by its highest-scoring match.
cat $SAMPLE_TYPE.single_exon_overlap.bed \
    $SAMPLE_TYPE.terminal_exon_overlap.bed \
    $SAMPLE_TYPE.upstream.bed \
    $SAMPLE_TYPE.internal_exon_overlap.bed \
    $SAMPLE_TYPE.intronic.bed \
    $SAMPLE_TYPE.antisense.bed \
    $SAMPLE_TYPE.anti_intronic.bed > $SAMPLE_TYPE.gene.bed
awk '{ printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t.\t.\t0\n"}' $SAMPLE_TYPE.peaks.bed >> $SAMPLE_TYPE.gene.bed
# Add a +-1 buffer to each feature to ensure that untemplated upstream nucleotides are included in the range
bedtools sort -i $SAMPLE_TYPE.gene.bed | awk '{ printf $1"\t"$2"\t"$3+1"\t"$4"\t"$5+1"\t"$6"\t"$7"\t"$8"\t"$9"\n"}' > $SAMPLE_TYPE.gene.sorted.bed

python $python_dir/bed_deduplicate.py -F 3 --select highscore --scoreline 8 $SAMPLE_TYPE.gene.sorted.bed | bedtools sort > $SAMPLE_TYPE.all.sorted.bed
python $python_dir/bed_deduplicate.py -F 3 --select highscore --scoreline 8 $SAMPLE_TYPE.all.sorted.bed \
    | sed 's/\t7$/\tS/' \
    | sed 's/\t6$/\tP/' \
    | sed 's/\t5$/\tD/' \
    | sed 's/\t4$/\tU/' \
    | sed 's/\t3$/\tI/' \
    | sed 's/\t2$/\tA/' \
    | sed 's/\t1$/\tIA/' \
    | sed 's/\t0$/\tN/' \
    > $SAMPLE_TYPE.all_features.unsorted.bed

bedtools sort -i $SAMPLE_TYPE.all_features.unsorted.bed | awk '{$4 = "5P."NR; print}' | sed 's/ /\t/g' > $SAMPLE_TYPE.all_features.bed

rm -f $SAMPLE_TYPE.single_exon_overlap.bed $SAMPLE_TYPE.gene.bed $SAMPLE_TYPE.gene.sorted.bed \
      $SAMPLE_TYPE.all.sorted.bed $SAMPLE_TYPE.all_features.unsorted.bed \
      $SAMPLE_TYPE.upstream.bed $SAMPLE_TYPE.terminal_exon_overlap.bed $SAMPLE_TYPE.is_upstream_of.bed \
      $SAMPLE_TYPE.internal_exon_overlap.bed $SAMPLE_TYPE.nearest_downstream.bed $SAMPLE_TYPE.residuals.bed \
      $SAMPLE_TYPE.nearest_upstream.bed $SAMPLE_TYPE.intronic.bed $SAMPLE_TYPE.all_sense.bed $SAMPLE_TYPE.antisense.bed \
      $SAMPLE_TYPE.anti_upstream_of.bed $SAMPLE_TYPE.anti_intronic.bed $SAMPLE_TYPE.anti_downstream_of.bed

########################################################


echo "Splitting capped and noncapped features..."

# Expected input format for bed_uug_filter.py:
# chromosome    start   stop    name    distance    strand  peak_location    closest_gene_ID

python $python_dir/bed_uug_filter.py \
    -C $SAMPLE_TYPE.capped.bed \
    -U $SAMPLE_TYPE.noncapped.bed \
    --min_uug $UUG \
    $SAMPLE_TYPE.all_features.bed \
    $SAMPLE_TYPE.plus.bedgraph \
    $SAMPLE_TYPE.minus.bedgraph \
    $SAMPLE_TYPE.uG.plus.bedgraph \
    $SAMPLE_TYPE.uG.minus.bedgraph

bedtools sort -i $SAMPLE_TYPE.capped.bed | awk '{ $4="5P.capped."NR; print }' | sed 's/ /\t/g' > $SAMPLE_TYPE.capped.sorted.bed
bedtools sort -i $SAMPLE_TYPE.noncapped.bed | awk '{ $4="5P.noncapped."NR; print }' | sed 's/ /\t/g' > $SAMPLE_TYPE.noncapped.sorted.bed
cat $SAMPLE_TYPE.capped.sorted.bed $SAMPLE_TYPE.noncapped.sorted.bed | bedtools sort > $SAMPLE_TYPE.5P_features.bed

rm -f $SAMPLE_TYPE.capped.sorted.bed $SAMPLE_TYPE.noncapped.sorted.bed \
      $SAMPLE_TYPE.all_features.bed $SAMPLE_TYPE.peaks.bed $SAMPLE_TYPE.only_peaks.bed \
      $SAMPLE_TYPE.rep.bed

echo "Cap masking bedgraph files..."

capped_bedgraphs=()
noncapped_bedgraphs=()
names=()

for s in ${samples[@]}
do
    python $mask \
    -P $endgraph_dir/$s/5P.plus_mask.bedgraph \
    -M $endgraph_dir/$s/5P.minus_mask.bedgraph \
    -PO $data_folder/$s.capmask.plus.bedgraph \
    -MO $data_folder/$s.capmask.minus.bedgraph \
    -I $SAMPLE_TYPE.capped.bed \
    -L $length_table
    
    names+=( $s.plus $s.minus )
    capped_bedgraphs+=( $endgraph_dir/$s/5P.plus_mask.bedgraph $endgraph_dir/$s/5P.minus_mask.bedgraph )
    noncapped_bedgraphs+=( $data_folder/$s.capmask.plus.bedgraph $data_folder/$s.capmask.minus.bedgraph )
done

echo "Cap masking merged bedgraph..."
python $mask \
-P $SAMPLE_TYPE.plus.bedgraph \
-M $SAMPLE_TYPE.minus.bedgraph \
-PO $SAMPLE_TYPE.capmask.plus.bedgraph \
-MO $SAMPLE_TYPE.capmask.minus.bedgraph \
-I $SAMPLE_TYPE.capped.bed \
-L $length_table

###########################
# CALCULATE READ COVERAGE #
###########################

# Subset 5' end features for only those that are
# (1) capped, and
# (2) within 500nt of an existing gene annotation

bedtools subtract -a $resource_dir/class.exons_by_gene.bed -b $SAMPLE_TYPE.capped.bed -s > $SAMPLE_TYPE.exons_noncapped.bed
grep -P '\t[PIUD]\t' $SAMPLE_TYPE.capped.bed |\
    awk '{ printf $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6"\n" }' > $SAMPLE_TYPE.exons_capped.bed


echo "Calculating gene-level capped and noncapped read coverage..."
python $python_dir/bed_feature_coverage.py \
    -L $length_table -F $SAMPLE_TYPE.exons_capped.bed \
    -I ${plus_list[@]} ${minus_list[@]} \
    -N ${plus_names[@]} ${minus_names[@]} \
    -O $SAMPLE_TYPE.capped.counts.tsv

python $python_dir/bed_feature_coverage.py \
    -L $length_table -F $SAMPLE_TYPE.exons_noncapped.bed \
    -I ${plus_list[@]} ${minus_list[@]} \
    -N ${plus_names[@]} ${minus_names[@]} \
    -O $SAMPLE_TYPE.noncapped.counts.tsv

rm -f $SAMPLE_TYPE.exons*capped.bed
for s in ${samples[@]}
do
    rm $s.*us.bedgraph
done



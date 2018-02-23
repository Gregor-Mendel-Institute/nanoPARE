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
endclass_dir=$results_dir/EndClass
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
WEIGHTED="true"
if [ $WEIGHTED == "true" ]
then
    ADJUSTMENT=( "R" )
else
    ADJUSTMENT=( "R" "W" )
fi

#####################
# ENVIRONMENT SETUP #
#####################

length_table=$resource_dir/length.table
mask=$python_dir/bedgraph_mask.py
coverage=$python_dir/bed_feature_coverage.py
samples=($sample_list)
data_folder=$output_folder/$sample_name

echo "Samples: ${samples[@]}"
echo "Sample type: $sample_name"

mkdir -p $data_folder
cd $data_folder

echo "Setting up exon reference..."
grep -P "exon\t" $annotation_gff | sed 's/\t[^\t]*transcript_id=\([^;]*\);.*/\t\1\t/' | sed 's/Parent=\([^;]*\).*/\1/' > exons_by_transcript.gff
python $python_dir/bed_deduplicate.py -F 8 -S upstream --startline 3 --endline 4 --strandline 6 exons_by_transcript.gff > terminal_exons_by_transcript.gff
sed 's/\.[0-9]$//' terminal_exons_by_transcript.gff > terminal_exons_by_gene.gff
grep -P "exon\t" $annotation_gff | sed 's/gene_id=\([^;]*\);.*/\1\t/' | sed 's/Parent=\([^;\.]*\).*/\1/' > exons_by_gene.gff

awk '{printf $1"\t"$4-1"\t"$5"\t"$9"\t"$6"\t"$7"\t"$8"\n"}' exons_by_gene.gff | bedtools sort > exons.bed
awk '{printf $1"\t"$4-1"\t"$5"\t"$9"\t"$6"\t"$7"\t"$8"\n"}' terminal_exons_by_gene.gff | bedtools sort > terminal_exons.bed

##################
# MERGE FEATURES #
##################

for A in ${ADJUSTMENT[@]}
do
    echo "Merging feature files..."
    rm -f $sample_name."$A".all.bed
    touch $sample_name."$A".all.bed
    plus_list=()
    minus_list=()
    plus_list_uug=()
    minus_list_uug=()
    for s in ${samples[@]}
    do
        sed "s/\t5P\./\t"$s"\.5P\./" \
            $endgraph_dir/$s/$s."$A".end_features.bed \
            >> $sample_name."$A".all.bed
        
        plus_list+=($endgraph_dir/$s/5P."$A"_plus_mask.bedgraph)
        minus_list+=($endgraph_dir/$s/5P."$A"_minus_mask.bedgraph)
        plus_list_uug+=($endgraph_dir/$s/uuG_plus_mask.bedgraph)
        minus_list_uug+=($endgraph_dir/$s/uuG_minus_mask.bedgraph)
    done
    
    echo "Bedgraph files +: "${plus_list[@]}
    echo "Bedgraph files -: "${minus_list[@]}
    python $python_dir/bedgraph_combine.py \
        -i ${plus_list[@]} \
        -o $sample_name."$A".plus.bedgraph
    
    python $python_dir/bedgraph_combine.py \
        -i ${minus_list[@]} \
        -o $sample_name."$A".minus.bedgraph
    
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
        -i $sample_name."$A".all.bed \
        > $sample_name."$A".merged.bed
    
    # Keep only replicable features (present in more than one library)
    grep -v -P '\t1\t' $sample_name."$A".merged.bed \
        > $sample_name."$A".rep.bed
    
    rm $sample_name."$A".all.bed $sample_name."$A".merged.bed
    
    bedfile=$sample_name."$A".rep.bed
    
    # Adding a merged peak position to each replicable feature, keeping both the closest gene ID and distance
    python $python_dir/bed_find_peaks.py \
        -I $bedfile \
        -O $sample_name."$A".peaks.bed \
        -BP $sample_name."$A".plus.bedgraph \
        -BM $sample_name."$A".minus.bedgraph \
        -L $length_table \
        -V pass
    
    bedtools sort -i $sample_name."$A".peaks.bed | awk '{printf $1"\t"$2"\t"$3"\t5P."NR"\t"$5"\t"$6"\t"$7"\n"}' > $sample_name."$A".rep_with_peaks.bed
    awk '{ printf $1"\t"$2+$7"\t"$2+$7+1"\t"$4"\t"$5"\t"$6"\t"$2"\t"$3"\t"$7"\n" }' $sample_name."$A".rep_with_peaks.bed > $sample_name."$A".rep_only_peaks.bed
    rm -f $sample_name."$A".peaks.bed \
          $sample_name."$A".rep_with_peaks.bed \
          $sample_name."$A".peaks.bed \
          
    
    
    # Heierarchically label each feature based on its relationship to reference:
    # PT (4 pts): peak position overlaps an annotated terminal exon
    # PI (3 pts): peak position overlaps an annotated internal exon
    # FT (2 pts): >=1nt of the feature overlaps an annotated terminal exon
    # FI (1 pt): >=1nt of the feature overlaps an annotated internal exon
    # O (0 pts): other peak, no overlap to any annotated genes
    
    bedfile=$sample_name."$A".rep_only_peaks.bed
    bedtools closest \
        -s -id \
        -D b \
        -t first \
        -a $bedfile \
        -b terminal_exons_by_gene.gff \
        > $sample_name.PT.tmp.bed
    
    grep -P '\t0$' $sample_name.PT.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t4\n" }' > $sample_name.PT.bed
    rm $sample_name.PT.tmp.bed
    
    bedtools closest \
        -s -id \
        -D b \
        -t first \
        -a $bedfile \
        -b exons_by_gene.gff \
        > $sample_name.PI.tmp.bed
    
    grep -P '\t0$' $sample_name.PI.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t3\n" }' > $sample_name.PI.bed
    rm $sample_name.PI.tmp.bed
    
    bedfile=$sample_name.rep_with_peaks.bed
    bedtools closest \
        -s -id \
        -D b \
        -t first \
        -a $bedfile \
        -b terminal_exons_by_gene.gff \
        > $sample_name.FT.tmp.bed
    
    grep -P '\t0$' $sample_name.FT.tmp.bed | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6"\t"$16"\t"$17"\t2\n"}' > $sample_name.FT.bed
    rm $sample_name.FT.tmp.bed
    
    bedtools closest \
        -s -id \
        -D b \
        -t first \
        -a $bedfile \
        -b exons_by_gene.gff \
        > $sample_name.FI.tmp.bed
    
    grep -P '\t0$' $sample_name.FI.tmp.bed | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6"\t"$16"\t"$17"\t1\n"}' > $sample_name.FI.bed
    grep -v -P '\t0$' $sample_name.FI.tmp.bed | awk '{printf $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6"\t"$16"\t"$17"\t0\n"}' > $sample_name.O.bed
    rm $sample_name.FI.tmp.bed
    
    # Collapse the decision tree into a single BED file, with each feature
    # represented by its highest-scoring match.
    cat $sample_name.PT.bed $sample_name.PI.bed $sample_name.FT.bed $sample_name.FI.bed $sample_name.O.bed | bedtools sort > $sample_name.gene.bed
    python $python_dir/bed_deduplicate.py -F 3 --select highscore --scoreline 8 $sample_name.gene.bed \
        | sed 's/\t0\t4$/\tPT/' \
        | sed 's/\t0\t3$/\tPI/' \
        | sed 's/\t0\t2$/\tFT/' \
        | sed 's/\t0\t1$/\tFI/' \
        | sed 's/\t0$//' \
        | awk '{ printf $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$6"\t"$5"\t"$7"\n" }' \
        > $sample_name."$A".all_features.bed
    
    rm -f $sample_name.PT.bed \
          $sample_name.PI.bed \
          $sample_name.FT.bed \
          $sample_name.FI.bed \
          $sample_name.O.bed \
          $sample_name.gene.bed \
          $sample_name."$A".rep_only_peaks.bed
    
    echo "Splitting capped and noncapped features..."
    
    # Expected input format for bed_uug_filter.py:
    # chromosome    start   stop    name    distance    strand  peak_location    closest_gene_ID
    
    python $python_dir/bed_uug_filter.py \
        -C $sample_name."$A".capped.bed \
        -U $sample_name."$A".noncapped.bed \
        --min_uug 0.1 \
        $sample_name."$A".all_features.bed \
        $sample_name."$A".plus.bedgraph \
        $sample_name."$A".minus.bedgraph \
        $sample_name.uuG.plus.bedgraph \
        $sample_name.uuG.minus.bedgraph
    
    echo "Cap masking bedgraph files..."
    
    capped_bedgraphs=()
    noncapped_bedgraphs=()
    names=()
    
    for s in ${samples[@]}
    do
        python $mask \
        -P $endgraph_dir/$s/5P."$A"_plus_mask.bedgraph \
        -M $endgraph_dir/$s/5P."$A"_minus_mask.bedgraph \
        -PO $data_folder/$s."$A".capmask.plus.bedgraph \
        -MO $data_folder/$s."$A".capmask.minus.bedgraph \
        -I $sample_name."$A".capped.bed \
        -L $length_table
        
        names+=( $s.plus $s.minus )
        capped_bedgraphs+=( $endgraph_dir/$s/5P."$A"_plus_mask.bedgraph $endgraph_dir/$s/5P."$A"_minus_mask.bedgraph )
        noncapped_bedgraphs+=( $data_folder/$s.capmask.plus.bedgraph $data_folder/$s.capmask.minus.bedgraph )
    done
    
    echo "Cap masking merged bedgraph..."
    python $mask \
    -P $sample_name."$A".plus.bedgraph \
    -M $sample_name."$A".minus.bedgraph \
    -PO $sample_name."$A".capmask.plus.bedgraph \
    -MO $sample_name."$A".capmask.minus.bedgraph \
    -I $sample_name."$A".capped.bed \
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
        $sample_name."$A".capped.bed \
        > $sample_name."$A".overlapping.capped.bed
    
    cat $sample_name."$A".overlapping.capped.bed exons.bed > "$sample_name"."$A"_gene_features.bed
    rm $sample_name."$A".overlapping.capped.bed
    
    echo "Calculating gene-level total and uncapped read coverage..."
    python $coverage \
        -F $sample_name."$A".capped.bed \
        -I ${capped_bedgraphs[@]} \
        -N ${names[@]} \
        -O "$sample_name"."$A"_capped_coverage.tsv \
        -L $length_table 
    
    python $coverage \
        -F $sample_name."$A".noncapped.bed \
        -I ${capped_bedgraphs[@]} \
        -N ${names[@]} \
        -O "$sample_name"."$A"_noncapped_coverage.tsv \
        -L $length_table
    
    python $coverage \
        -F "$sample_name"."$A"_gene_features.bed \
        -I ${capped_bedgraphs[@]} \
        -N ${names[@]} \
        -O "$sample_name"."$A"_total_coverage.tsv \
        -L $length_table
    
    python $coverage \
        -F "$sample_name"."$A"_gene_features.bed \
        -I ${noncapped_bedgraphs[@]} \
        -N ${names[@]} \
        -O "$sample_name"."$A"_capmasked_coverage.tsv \
        -L $length_table
done


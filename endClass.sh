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
if [[ $WEIGHTED == "true" ]]
then
    ADJUSTMENT=( "R" "W" )
else
    ADJUSTMENT=( "R" )
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

if [ $SETUP == "true" ]
then
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
    
    echo "Setup complete. Annotation files stored in resource directory."
    exit 0
fi

##################
# MERGE FEATURES #
##################

for A in ${ADJUSTMENT[@]}
do
    echo "Merging feature files..."
    rm -f $sample_name."$A".all.bed
    touch $sample_name."$A".all.bed
    plus_list=()
    plus_names=()
    minus_list=()
    minus_names=()
    plus_list_uug=()
    minus_list_uug=()
    for s in ${samples[@]}
    do
        sed "s/\t5P\./\t"$s"\.5P\./" \
            $endgraph_dir/$s/$s."$A".end_features.bed \
            >> $sample_name."$A".all.bed
        
        plus_list+=($endgraph_dir/$s/5P."$A"_plus_mask.bedgraph)
        plus_names+=($s.plus)
        minus_list+=($endgraph_dir/$s/5P."$A"_minus_mask.bedgraph)
        minus_names+=($s.minus)
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
        -L $length_table --trim 1 \
        -V pass
    
    # Classify features against the annotated list of genes
    # 1: Generate a bed file of 1-nt features representing only the peaks of each 5P feature
    awk '{ printf $1"\t"$2+$5"\t"$2+$5+1"\t"$4"\t"$5"\t"$6"\t"$2"\t"$3"\t"$5"\n" }'\
        $sample_name."$A".peaks.bed | bedtools sort > $sample_name.$A.only_peaks.bed
    bedfile=$sample_name.$A.only_peaks.bed
    
    # 2: Find the nearest overlapping/downstream exon element in the single-exon transcripts
    bedtools closest -s -id -D b -t first \
        -a $bedfile \
        -b $resource_dir/single_exon_genes.bed > $sample_name.single.tmp.bed
    grep -P '\t0$' $sample_name.single.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$13"\t"$14"\t7\n" }' > $sample_name.single_exon_overlap.bed
    rm $sample_name.single.tmp.bed
    
    # 3: Find the nearest overlapping/downstream exon element in the annotation set
    bedtools closest -s -id -D b -t first \
        -a $bedfile \
        -b $resource_dir/terminal_exons_by_gene.gff > $sample_name.terminal.tmp.bed
    grep -P '\t0$' $sample_name.terminal.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t6\n" }' > $sample_name.terminal_exon_overlap.bed
    grep -v -P '\t0$' $sample_name.terminal.tmp.bed | awk '{ if ( $19 >= -500 ) printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t4\n" }' > $sample_name.upstream.bed
    rm $sample_name.terminal.tmp.bed

    bedtools closest -s -id -D b -t first \
        -a $bedfile \
        -b $resource_dir/exons_by_gene.gff > $sample_name.nearest_downstream.tmp.bed
    grep -P '\t0$' $sample_name.nearest_downstream.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t5\n" }' > $sample_name.exon_overlap.bed
    bedtools intersect -v -wa -a $sample_name.exon_overlap.bed -b $sample_name.terminal_exon_overlap.bed > $sample_name.internal_exon_overlap.bed
    grep -v -P '\t0$' $sample_name.nearest_downstream.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t0\n" }' > $sample_name.is_upstream_of.bed
    bedtools intersect -v -wa -a $sample_name.is_upstream_of.bed -b $sample_name.upstream.bed > $sample_name.nearest_downstream.bed
    rm $sample_name.nearest_downstream.tmp.bed $sample_name.exon_overlap.bed

    bedtools closest -s -iu -D b -t first \
        -a $bedfile \
        -b $resource_dir/exons_by_gene.gff > $sample_name.nearest_upstream.tmp.bed
    grep -v -P '\t0$' $sample_name.nearest_upstream.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t0\n" }' > $sample_name.is_downstream_of.bed
    bedtools intersect -v -wa -a $sample_name.is_downstream_of.bed -b $sample_name.upstream.bed > $sample_name.nearest_upstream.bed
    rm $sample_name.nearest_upstream.tmp.bed $sample_name.is_downstream_of.bed

    bedtools closest -a $sample_name.nearest_upstream.bed -b $sample_name.nearest_downstream.bed \
        | awk '{ if ( $18 == 0 && $7 == $16 ) printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$17"\t3\n" }' \
        > $sample_name.intronic.bed

    cat $sample_name.terminal_exon_overlap.bed $sample_name.upstream.bed $sample_name.internal_exon_overlap.bed $sample_name.intronic.bed \
        | bedtools sort > $sample_name.all_sense.bed
    bedtools intersect -v -s -wa -a $bedfile -b $sample_name.all_sense.bed > $sample_name.residuals.bed

    bedtools closest -S -id -D b -t first \
        -a $sample_name.residuals.bed \
        -b $resource_dir/exons_by_gene.gff > $sample_name.antisense_a.tmp.bed
    grep -P '\t0$' $sample_name.antisense_a.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t2\n" }' > $sample_name.antisense.bed
    grep -v -P '\t0$' $sample_name.antisense_a.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t0\n" }' > $sample_name.anti_upstream_of.bed

    bedtools closest -S -iu -D b -t first \
        -a $sample_name.residuals.bed \
        -b $resource_dir/exons_by_gene.gff > $sample_name.antisense_d.tmp.bed
    grep -v -P '\t0$' $sample_name.antisense_d.tmp.bed | awk '{ printf $1"\t"$7"\t"$8"\t"$4"\t"$9"\t"$6"\t"$18"\t"$19"\t0\n" }' > $sample_name.anti_downstream_of.bed
    bedtools closest -a $sample_name.anti_upstream_of.bed -b $sample_name.anti_downstream_of.bed \
        | awk '{ if ( $18 == 0 && $7 == $16 ) printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$17"\t1\n" }' \
        > $sample_name.anti_intronic.bed

    rm $sample_name.antisense_a.tmp.bed $sample_name.antisense_d.tmp.bed

    # Collapse the decision tree into a single BED file, with each feature
    # represented by its highest-scoring match.
    cat $sample_name.single_exon_overlap.bed \
        $sample_name.terminal_exon_overlap.bed \
        $sample_name.upstream.bed \
        $sample_name.internal_exon_overlap.bed \
        $sample_name.intronic.bed \
        $sample_name.antisense.bed \
        $sample_name.anti_intronic.bed > $sample_name.gene.bed
    awk '{ printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t.\t.\t0\n"}' $sample_name."$A".peaks.bed >> $sample_name.gene.bed
    # Add a +-1 buffer to each feature to ensure that untemplated upstream nucleotides are included in the range
    bedtools sort -i $sample_name.gene.bed | awk '{ printf $1"\t"$2-1"\t"$3+1"\t"$4"\t"$5+1"\t"$6"\t"$7"\t"$8"\t"$9"\n"}' > $sample_name.gene.sorted.bed
    
    python $python_dir/bed_deduplicate.py -F 3 --select highscore --scoreline 8 $sample_name.gene.sorted.bed | bedtools sort > $sample_name.all.sorted.bed
    python $python_dir/bed_deduplicate.py -F 3 --select highscore --scoreline 8 $sample_name.all.sorted.bed \
        | sed 's/\t7$/\tS/' \
        | sed 's/\t6$/\tP/' \
        | sed 's/\t5$/\tD/' \
        | sed 's/\t4$/\tU/' \
        | sed 's/\t3$/\tI/' \
        | sed 's/\t2$/\tA/' \
        | sed 's/\t1$/\tIA/' \
        | sed 's/\t0$/\tN/' \
        > $sample_name."$A".all_features.unsorted.bed

    bedtools sort -i $sample_name."$A".all_features.unsorted.bed | awk '{$4 = "5P."NR; print}' | sed 's/ /\t/g' > $sample_name."$A".all_features.bed
    
    rm -f $sample_name.single_exon_overlap.bed $sample_name.gene.bed $sample_name.gene.sorted.bed $sample_name.all.sorted.bed $sample_name.all_features.unsorted.bed \
        $sample_name.upstream.bed $sample_name.terminal_exon_overlap.bed $sample_name.is_upstream_of.bed \
        $sample_name.internal_exon_overlap.bed $sample_name.nearest_downstream.bed $sample_name.residuals.bed \
        $sample_name.nearest_upstream.bed $sample_name.intronic.bed $sample_name.all_sense.bed $sample_name.antisense.bed \
        $sample_name.anti_upstream_of.bed $sample_name.anti_intronic.bed $sample_name.anti_downstream_of.bed
    
    ########################################################
    
    
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
    
    bedtools sort -i $sample_name."$A".capped.bed | awk '{ $4="5P.capped."NR; print }' | sed 's/ /\t/g' > $sample_name."$A".capped.sorted.bed
    bedtools sort -i $sample_name."$A".noncapped.bed | awk '{ $4="5P.noncapped."NR; print }' | sed 's/ /\t/g' > $sample_name."$A".noncapped.sorted.bed
    cat $sample_name."$A".capped.sorted.bed $sample_name."$A".noncapped.sorted.bed | bedtools sort > $sample_name."$A".5P.bed
    
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
        noncapped_bedgraphs+=( $data_folder/$s."$A".capmask.plus.bedgraph $data_folder/$s."$A".capmask.minus.bedgraph )
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
    
    bedtools subtract -a $resource_dir/exons_by_gene.bed -b $sample_name."$A".capped.bed -s > $sample_name."$A".exons_noncapped.bed
    grep -P '\t[PIUD]\t' $sample_name."$A".capped.bed |\
        awk '{ printf $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6"\n" }' > $sample_name."$A".exons_capped.bed
    
    
    echo "Calculating gene-level capped and noncapped read coverage..."
    python $python_dir/bed_feature_coverage.py \
        -L $length_table -F $sample_name."$A".exons_capped.bed \
        -I ${plus_list[@]} ${minus_list[@]} \
        -N ${plus_names[@]} ${minus_names[@]} \
        -O $sample_name."$A".capped.counts.tsv
    
    python $python_dir/bed_feature_coverage.py \
        -L $length_table -F $sample_name."$A".exons_noncapped.bed \
        -I ${plus_list[@]} ${minus_list[@]} \
        -N ${plus_names[@]} ${minus_names[@]} \
        -O $sample_name."$A".noncapped.counts.tsv
        
done


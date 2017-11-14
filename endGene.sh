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
    reference_table=$resource_dir/reference_table_EndGene.txt
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
    output_folder=$results_dir/EndGene
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
MASK_NAME=${input_array[3]} # name of the sample type whose capped features should be used

#####################
# ENVIRONMENT SETUP #
#####################

length_table=$resource_dir/length.table
tso_mask_file=$resource_dir/transcript_TSO_mask.bed
mask=$python_dir/bedgraph_mask.py
coverage=$python_dir/bed_feature_coverage.py
endmap_folder=$results_dir/EndMap
endmask_folder=$results_dir/EndMask

endgraph_folder=$temp_dir
output_folder=$results_dir/EndGene
data_folder=$output_folder/$SAMPLE_NAME
capped_features=$endmask_folder/$MASK_NAME/$MASK_NAME.capped.bed

rm -R $data_folder
mkdir -p $data_folder
cd $data_folder

samples=($sample_list)
echo "Samples: ${samples[@]}"
echo "Sample type: $SAMPLE_NAME"

###########################
# CALCULATE READ COVERAGE #
###########################

echo "Determining dominant transcript for $SAMPLE_NAME..."
python $python_dir/bedgraph_genome_to_transcripts.py \
    --subset $annotation_subset \
    --output $SAMPLE_NAME.all.transcript.bedgraph \
    $endmask_folder/$MASK_NAME/$MASK_NAME.plus.bedgraph \
    $endmask_folder/$MASK_NAME/$MASK_NAME.minus.bedgraph \
    $annotation_gff \
    $genome_fasta

python $python_dir/bedgraph_dominant_transcripts.py \
    -I $SAMPLE_NAME.all.transcript.bedgraph \
    -F $genome_fasta \
    -O $SAMPLE_NAME.dom.transcript.bedgraph \
    -L $SAMPLE_NAME.dominant_transcript_lengths.tsv

rm $SAMPLE_NAME.all.transcript.bedgraph $SAMPLE_NAME.dom.transcript.bedgraph

echo "Making gene-level bed files..."
bedgraph_list=()
for s in ${samples[@]}
do
    # Perform cap masking on all 5' end bedgraphs (genome level)
    python $mask \
        -P $endmap_folder/$s/"$s"_plus.5p.bedgraph \
        -M $endmap_folder/$s/"$s"_minus.5p.bedgraph \
        -PO $data_folder/$s.capmask.plus.bedgraph \
        -MO $data_folder/$s.capmask.minus.bedgraph \
        -I $capped_features \
        -L $length_table
    
    # Generate dominant transcript-level cap masked bedgraphs
    python $python_dir/bedgraph_genome_to_transcripts.py \
        --subset $SAMPLE_NAME.dominant_transcript_lengths.tsv \
        --output $s.transcript.capmasked.bedgraph \
        $data_folder/$s.capmask.plus.bedgraph \
        $data_folder/$s.capmask.minus.bedgraph \
        $annotation_gff \
        $genome_fasta
    
    # Additionally mask TSO sequences in transcript-level bedgraphs
    python $mask \
        -P --output $s.transcript.capmasked.bedgraph \
        -PO $data_folder/$s.transcript.bedgraph \
        -U $tso_mask_file \
        -L $SAMPLE_NAME.dominant_transcript_lengths.tsv
    
    # Remove intermediate files
    rm $s.transcript.capmasked.bedgraph
    bedgraph_list+=( $s.transcript.bedgraph )
done

bedtools unionbedg -i ${bedgraph_list[@]} > $SAMPLE_NAME.merge.transcript.bedgraph

colnum=$(head -n 1 $SAMPLE_NAME.merge.transcript.bedgraph | awk '{print NF}')
i=4
sumstring='$4'
while [ "$i" -lt "$colnum" ]
do
  i=$(($i + 1))
  sumstring+='+$'$i
done

awk '{printf $1"\t"$2"\t"$3"\t"'"$sumstring"'"\n"}'  $SAMPLE_NAME.merge.transcript.bedgraph > $SAMPLE_NAME.transcript.bedgraph
rm $SAMPLE_NAME.merge.transcript.bedgraph

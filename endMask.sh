#!/bin/bash

function usage() {
cat <<HELP

### EndMask ###
Usage: ./endMask.sh [options] -T|--type <sample_type>

Masks reads belonging to capped features in the genome and prepares transcript bedgraph files for EndCut.
 1) Identifies the transcript isoform with the most reads from each gene
 2) Sets values for all 5P reads within a capped feature to 0
 3) Writes a transcript-indexed bedgraph file of cap-masked reads with all dominant transcripts

Expects EndMap, EndGraph, and EndClass to be run before and their results to be in directories
/results/EndMap/sample_name
/results/EndGraph/sample_name
/results/EndClass/sample_type
respectively, for all sample_names of the chosen sample_type that appear in reference.table

Optional arguments:
-R | --reference     Reference table (default: resources/reference.table)
-G | --genome        Genome FASTA file (default: resources/genome.fasta)
-A | --annotation    Transcript GFF file (default: resources/annotation.gff)
--lmod               Load required modules with Lmod (default: false)
--cpus               Number of cores available for multithreaded programs (default: 1)
--mask               Alternative sample_type to use for cap masking (default: sample_type)

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
r_dir=$root_dir/scripts/r_scripts
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
    usage
    exit 1
else
    . $bash_dir/read_cmdline.sh
fi

# Set defaults for variables if they are not already in the environment
if [ -z "$SAMPLE_TYPE" ]
then
    echo "ERROR: Please input a sample type to process from reference.table (-T|--type)"
    exit 1
fi
if [ -z "$MASK_NAME" ]
then
    MASK_NAME=$SAMPLE_TYPE
fi
if [ -z "$JOB_NUMBER" ]
then
    JOB_NUMBER=${PBS_ARRAY_INDEX} # Imports the PBS job array number if it exists. Can be overridden with the commandline argument -J $JOB_NUMBER
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
if [ -z "$transcript_fasta" ]
then
    transcript_fasta=$resource_dir/transcriptome.fasta
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

echo "Config settings:"
. $bash_dir/list_settings.sh

# Environment modules to load with Lmod (if option --lmod is passed)
REQUIRED_MODULES=( --bedtools )
. $bash_dir/load_modules.sh
echo " "

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

#####################
# ENVIRONMENT SETUP #
#####################

length_table=$resource_dir/length.table
mask=$python_dir/bedgraph_mask.py
coverage=$python_dir/bed_feature_coverage.py
data_folder=$output_folder/$SAMPLE_TYPE

rm -rf $data_folder
mkdir -p $data_folder
cd $data_folder

echo "Samples: ${samples[@]}"
echo "Sample type: $SAMPLE_TYPE"

###########################
# CALCULATE READ COVERAGE #
###########################

capped_features=$endclass_dir/$MASK_NAME/$MASK_NAME.capped.bed
echo "Determining dominant transcripts for $MASK_NAME..."
python $python_dir/bedgraph_genome_to_transcripts.py \
    --output $SAMPLE_NAME.all.transcript.bedgraph \
    $endclass_dir/$MASK_NAME/$MASK_NAME.plus.bedgraph \
    $endclass_dir/$MASK_NAME/$MASK_NAME.minus.bedgraph \
    $ANNOTATION_GFF \
    $GENOME_FASTA

python $python_dir/bedgraph_dominant_transcripts.py \
    -I $SAMPLE_NAME.all.transcript.bedgraph \
    -F $transcript_fasta \
    -O $SAMPLE_NAME.dom.transcript.bedgraph \
    -L $SAMPLE_NAME.dominant_transcript_lengths.tsv

rm $SAMPLE_NAME.all.transcript.bedgraph $SAMPLE_NAME.dom.transcript.bedgraph

echo "Making gene-level bed files..."
for s in ${samples[@]}
do
    # Perform cap masking on all 5' end bedgraphs (genome level)
    echo Capmasking: $s
    python $mask \
        -P $endmap_dir/$s/"$s".5P.plus.bedgraph \
        -M $endmap_dir/$s/"$s".5P.minus.bedgraph \
        -PO $data_folder/$s.capmask.plus.bedgraph \
        -MO $data_folder/$s.capmask.minus.bedgraph \
        -I $capped_features \
        -L $length_table
        
    # Generate dominant transcript-level cap masked bedgraphs
    echo Converting to transcripts: $s
    python $python_dir/bedgraph_genome_to_transcripts.py \
        --subset $SAMPLE_NAME.dominant_transcript_lengths.tsv \
        --output $s.transcript.capmasked.bedgraph \
        $data_folder/$s.capmask.plus.bedgraph \
        $data_folder/$s.capmask.minus.bedgraph \
        $ANNOTATION_GFF \
        $GENOME_FASTA
    
done
echo "EndMask complete! See results in results/EndMask/$SAMPLE_TYPE"

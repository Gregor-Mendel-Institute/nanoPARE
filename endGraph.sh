#!/bin/bash

function usage() {
cat <<HELP

### EndGraph ###
Usage: ./endGraph.sh [options] -N|--name <sample_name>

Identification of end features through subtractive kernel density estimation
 1) Determines a scaling factor to adjust read depths of 5P and BODY libraries
 2) Smooths signal by fitting a Laplace kernel to END - BODY read values
 3) Converts continuous regions of positive signal to features in a BED file

Expects EndMap to be run and the results for sample_name in the directory /results/EndMap/sample_name

Optional arguments:
-R | --reference     Reference table (default: resources/reference.table)
-G | --genome        Genome FASTA file (default: resources/genome.fasta)
-A | --annotation    Transcript GFF file (default: resources/annotation.gff)
--lmod               Load required modules with Lmod (default: false)
--cpus               Number of cores available for multithreaded programs (default: 1)
--rpm                Minimum RPM required to keep a feature (default: 0.5)
--kernel             Type of kernel to use for density estimation (default: laplace. options: gaussian, laplace)
--bandwidth          Bandwidth of kernel to use, in nucleotides (default: 15)
--fraglen            Mean fragment length of cDNA library, in nucleotides (default: 200)

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
genome_index_dir=$temp_dir/genome_index
ENDMAP_DIR=$results_dir/EndMap
ENDGRAPH_DIR=$results_dir/EndGraph

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
if [ -z "$SAMPLE_NAME" ]
then
    if [ -z "$SAMPLE_ARRAY" ]
    then
        echo "ERROR: Please input a sample name (-N|--name)"
        exit 1
    fi
    SAMPLES=( $(echo $SAMPLE_ARRAY | tr '!' ' ' ) )
    SAMPLE_NAME=${SAMPLES[$JOB_NUMBER]}
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
if [ -z "$FRAGLEN" ]
then
    FRAGLEN=200
fi
if [ -z "$BANDWIDTH" ]
then
    BANDWIDTH=15
fi
if [ -z "$RPM" ]
then
    RPM=0.5
fi
if [ -z "$KERNEL" ]
then
    KERNEL='laplace'
fi

echo "Config settings:"
. $bash_dir/list_settings.sh

# Environment modules to load with Lmod (if option --lmod is passed)
REQUIRED_MODULES=( --bedtools --python )
. $bash_dir/load_modules.sh
echo " "

sample_dir=$ENDGRAPH_DIR/$SAMPLE_NAME
mkdir -p $sample_dir

### STEP 1: BEDGRAPH ARTIFACT MASKING ###
echo "#########################################"
echo "### STEP 1: BEDGRAPH ARTIFACT MASKING ###"
echo "#########################################"
echo " "
# Generates a list of BED features indicating sites in the genome 
# where 5P or 3P reads are likely methodological artifacts:
#   TSO oligos perform strand invasion at internal sites complementary
#     to their 3' ends (Tang et al. 2012, Nucleic Acids Research)
#   oligo(dT) priming results in internal priming artifacts at templated 
#     A-rich or purine-rich regions of transcripts
# Customize what sequences are masked by editing the file 
# /resources/mask_sequences.table

#TODO: Update behavior for $library_type == 5P or 3P
library_type="5P"

# Copy background bedgraph files from results/EndMap
cat $ENDMAP_DIR/$SAMPLE_NAME/"$SAMPLE_NAME".BODY.plus.bedgraph > $sample_dir/BODY.plus.bedgraph
cat $ENDMAP_DIR/$SAMPLE_NAME/"$SAMPLE_NAME".BODY.minus.bedgraph > $sample_dir/BODY.minus.bedgraph

# Copy 5P/3P begraph files from results/EndMap
if [ $library_type == "5P" ]
then
    python $python_dir/bedgraph_mask.py \
        -P $ENDMAP_DIR/$SAMPLE_NAME/"$SAMPLE_NAME".5P.plus.bedgraph \
        -M $ENDMAP_DIR/$SAMPLE_NAME/"$SAMPLE_NAME".5P.minus.bedgraph \
        -PO $sample_dir/5P.plus_mask.bedgraph \
        -MO $sample_dir/5P.minus_mask.bedgraph \
        -U $resource_dir/5P_mask_up.bed \
        -L $resource_dir/length.table
    
    python $python_dir/bedgraph_mask.py \
        -P $ENDMAP_DIR/$SAMPLE_NAME/"$SAMPLE_NAME".uG.plus.bedgraph \
        -M $ENDMAP_DIR/$SAMPLE_NAME/"$SAMPLE_NAME".uG.minus.bedgraph \
        -PO $sample_dir/uG.plus_mask.bedgraph \
        -MO $sample_dir/uG.minus_mask.bedgraph \
        -U $resource_dir/5P_mask_up.bed \
        -L $resource_dir/length.table
fi
if [ $library_type == "3P" ]
then
    python $python_dir/bedgraph_mask.py \
        -P $sample_dir/"$SAMPLE_NAME".3P.plus.bedgraph \
        -M $sample_dir/"$SAMPLE_NAME".3P.minus.bedgraph \
        -PO $sample_dir/3P.plus_mask.bedgraph \
        -MO $sample_dir/3P.minus_mask.bedgraph \
        -D $resource_dir/3P_mask_down.bed \
        -L $resource_dir/length.table
fi

echo "Step 1 complete."

### STEP 2: PARAMETER ESTIMATION WITH SUBTRACTIVE BEDGRAPH FILES ###
echo "####################################################################"
echo "### STEP 2: PARAMETER ESTIMATION WITH SUBTRACTIVE BEDGRAPH FILES ###"
echo "####################################################################"
echo " "
# Performs a meta-analysis of annotated transcripts to determine:
#     Scaling factor to adjust for capture efficiency of 5P/3P/BODY reads

scale=1
SCALE_CAP=10
BANDWIDTH_CAP=30

flip_5p="true"
if [[ $library_type == "5P" ]]
then
    if [[ $flip_5p == "true" ]]
    then
        bg_plus=$sample_dir/BODY.minus.bedgraph
        bg_minus=$sample_dir/BODY.plus.bedgraph
    else
        bg_plus=$sample_dir/BODY.plus.bedgraph
        bg_minus=$sample_dir/BODY.minus.bedgraph    
    fi
else
    bg_plus=$sample_dir/BODY.plus.bedgraph
    bg_minus=$sample_dir/BODY.minus.bedgraph
fi

# Write quantification tables for scaling factor estimation
python $python_dir/gtf_quantify.py \
    -A $ANNOTATION_GFF \
    -P $sample_dir/"$library_type".plus_mask.bedgraph \
    -M $sample_dir/"$library_type".minus_mask.bedgraph \
    --norm reads length RPM TPM \
    --buffer 100 \
    --gene > $sample_dir/"$SAMPLE_NAME"."$library_type".quant.tsv

python $python_dir/gtf_quantify.py \
    -A $ANNOTATION_GFF \
    -P $bg_plus \
    -M $bg_minus \
    --norm reads length RPM TPM \
    --buffer 100 \
    --gene > $sample_dir/"$SAMPLE_NAME".bg.quant.tsv

# Perform scaling based on quantification tables above
python $python_dir/endgraph_ratio.py \
    $sample_dir/"$SAMPLE_NAME"."$library_type".quant.tsv \
    $sample_dir/"$SAMPLE_NAME".bg.quant.tsv \
    -F $FRAGLEN > $sample_dir/scaling_factor.txt

readnum_5P=$(cat $sample_dir/"$SAMPLE_NAME".5P.quant.tsv | awk '{ sum += $2 } END { print sum }')
readnum_bg=$(cat $sample_dir/"$SAMPLE_NAME".bg.quant.tsv | awk '{ sum += $2 } END { print sum }')

scale=$(cat $sample_dir/scaling_factor.txt)
echo "Gene-mapping 5P reads:   $readnum_5P"
echo "Gene-mapping BODY reads: $readnum_bg"
echo "Scaling factor: $scale"
echo " "

if [ $(echo "$scale > $SCALE_CAP" | bc -l) -eq 1 ]
then    
    scale=$SCALE_CAP
    echo "Hit scale cap of $SCALE_CAP"
fi

if [ $(echo "$BANDWIDTH > $BANDWIDTH_CAP" | bc -l) -eq 1 ]
then
    BANDWIDTH=$BANDWIDTH_CAP
    echo "Cannot exceed bandwidth cap of $BANDWIDTH_CAP"
fi

python $python_dir/bedgraph_combine.py \
    -i $sample_dir/"$library_type".plus_mask.bedgraph $bg_plus  \
    -s $scale -1 \
    -o $sample_dir/"$library_type".plus_subtract.bedgraph

python $python_dir/bedgraph_combine.py \
    -i $sample_dir/"$library_type".minus_mask.bedgraph $bg_minus \
    -s $scale -1 \
    -o $sample_dir/"$library_type".minus_subtract.bedgraph

### STEP 3: CONTINUOUS KERNEL DENSITY DISTRIBUTION ###
echo "#########################################################"
echo "### STEP 3: CONTINUOUS KERNEL DENSITY DISTRIBUTION ###"
echo "#########################################################"
echo " "
# Takes the two subtractive BEDGRAPH files generated in Step 2
# and smooths them using a continous kernel function.
# See the Python util "bedgraph_kernel_density.py" for full details.
# Defaults to a summed Laplace distribution smoothing

for strand in plus minus
do
    kernel_density_command="python \
        $python_dir/bedgraph_kernel_density.py \
        -B $sample_dir/"$library_type"."$strand"_subtract.bedgraph \
        -O $sample_dir/"$library_type"."$strand"_smooth.bedgraph \
        -L $resource_dir/length.table \
        -K $KERNEL \
        -H $BANDWIDTH \
        -S 3 \
        -D 3 \
        -P \
        -c $CPUS"
    echo $kernel_density_command
    eval $kernel_density_command
    if [ ! -f $sample_dir/"$library_type"."$strand"_smooth.bedgraph ]
    then
        echo "ERROR: Failed to generate "$library_type"."$strand"_smooth.bedgraph"
        exit 1
    fi
    feature_threshold_command="python \
    $python_dir/bedgraph_thresh_to_bed.py \
    -B $sample_dir/"$library_type"."$strand"_smooth.bedgraph \
    -O $sample_dir/"$library_type"."$strand"_features.bed \
    -L $resource_dir/length.table \
    -T 0 \
    -M 10 \
    -V sum \
    -S $strand"
    echo $feature_threshold_command
    eval $feature_threshold_command
    if [ ! -f $sample_dir/"$library_type"."$strand"_features.bed ]
    then
        echo "ERROR: Failed to generate "$library_type"."$strand"_features.bed"
        exit 1
    fi
done

# Merges all end features identified in PHASE 3.4 to a single BED file.
rm -f $sample_dir/end_features_temp.bed
touch $sample_dir/end_features_temp.bed
sed 's/thresh./'$library_type'.plus./' $sample_dir/"$library_type".plus_features.bed \
    >> $sample_dir/end_features_temp.bed
sed 's/thresh./'$library_type'.minus./' $sample_dir/"$library_type".minus_features.bed \
    >> $sample_dir/end_features_temp.bed

bedtools sort -i $sample_dir/end_features_temp.bed > $sample_dir/end_features.bed

python $python_dir/bed_find_peaks.py \
    -I $sample_dir/end_features.bed \
    -L $resource_dir/length.table \
    -O $sample_dir/$SAMPLE_NAME.rpm.features.bed \
    -BP $sample_dir/"$library_type".plus_mask.bedgraph \
    -BM $sample_dir/"$library_type".minus_mask.bedgraph \
    -V rpm

awk -F'[\t]' -v rpm="$RPM" '{if ($7 >= rpm){ print }}' $sample_dir/$SAMPLE_NAME.rpm.features.bed > $sample_dir/$SAMPLE_NAME.end_features.bed

rm $sample_dir/end_features_temp.bed \
    $sample_dir/end_features.bed \
    $sample_dir/$SAMPLE_NAME.rpm.features.bed \
    $sample_dir/"$library_type".plus_features.bed \
    $sample_dir/"$library_type".minus_features.bed \
    $sample_dir/"$library_type".plus_subtract.bedgraph \
    $sample_dir/"$library_type".minus_subtract.bedgraph

echo "Step 3 complete."
echo "EndGraph complete!"


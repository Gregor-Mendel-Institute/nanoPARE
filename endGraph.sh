#!/bin/bash

# EndGraph
# Identification of end features through subtractive kernel density estimation
# 1) Determines a scaling factor to adjust read depths of 5P/3P and BODY libraries
# 2) Signal smoothing by fitting a kernel to END - BODY read values
# 3) Convert continuous regions of positive signal to features in a BED file
# Expects EndMap to be run and the results in the directory /results/EndMap/

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

sample_dir=$temp_dir/$sample_name
mkdir -p $sample_dir

if [ $SETUP == "true" ]
then
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
    echo "Setup complete."
    exit 0
fi

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

cp $temp_dir/$body_name/"$body_name"_plus.bedgraph $sample_dir/BODY_plus.bedgraph
cp $temp_dir/$body_name/"$body_name"_minus.bedgraph $sample_dir/BODY_minus.bedgraph

if [ $library_type == "5P" ]
then
    python $python_dir/bedgraph_mask.py \
        -P $sample_dir/"$sample_name"_plus.bedgraph \
        -M $sample_dir/"$sample_name"_minus.bedgraph \
        -PO $sample_dir/5P_plus_mask.bedgraph \
        -MO $sample_dir/5P_minus_mask.bedgraph \
        -U $resource_dir/5P_mask_up.bed \
        -L $resource_dir/length.table
    
    python $python_dir/bedgraph_mask.py \
        -P $sample_dir/"$sample_name"_plus.uG.bedgraph \
        -M $sample_dir/"$sample_name"_minus.uG.bedgraph \
        -PO $sample_dir/uuG_plus_mask.bedgraph \
        -MO $sample_dir/uuG_minus_mask.bedgraph \
        -U $resource_dir/5P_mask_up.bed \
        -L $resource_dir/length.table
fi
if [ $library_type == "3P" ]
then
    python $python_dir/bedgraph_mask.py \
        -P $sample_dir/"$sample_name"_plus.bedgraph \
        -M $sample_dir/"$sample_name"_minus.bedgraph \
        -PO $sample_dir/3P_plus_mask.bedgraph \
        -MO $sample_dir/3P_minus_mask.bedgraph \
        -D $resource_dir/3P_mask_down.bed \
        -L $resource_dir/length.table
fi

echo "Step 1 complete."

### STEP 2: PARAMETER ESTIMATION WITH SUBTRACTIVE BEDGRAPH FILES ###
echo "####################################################################"
echo "### STEP 2: PARAMETER ESTIMATION WITH SUBTRACTIVE BEDGRAPH FILES ###"
echo "####################################################################"
echo " "
# Performs a meta-analysis of annotated transcripts to determine 2 parameters:
#     Optimal bandwidth for kernel density estimation
#     Scaling factor to adjust for capture efficiency of 5P/3P/BODY reads
#
# Generates four scaled subtractive bedgraph files:
#     5Pplus  - BODYminus
#     5Pminus - BODYplus
#     3Pplus  - BODYplus
#     3Pminus - BODYminus

scale=1
SCALE_CAP=10
BANDWIDTH_CAP=30

flip_5p="true"

if [[ $library_type == "5P" ]]
then
    if [[ $flip_5p == "true" ]]
    then
        bg_plus=$sample_dir/BODY_minus.bedgraph
        bg_minus=$sample_dir/BODY_plus.bedgraph
    else
        bg_plus=$sample_dir/BODY_plus.bedgraph
        bg_minus=$sample_dir/BODY_minus.bedgraph    
    fi
else
    bg_plus=$sample_dir/BODY_plus.bedgraph
    bg_minus=$sample_dir/BODY_minus.bedgraph
fi


python $python_dir/bedgraph_rescale.py \
    -P $sample_dir/"$library_type"_plus_mask.bedgraph \
    -N $bg_plus \
    -A $annotation_gff \
    -S + \
    --position $library_type \
    --align \
    > $sample_dir/"$library_type"_scale_plus.txt

python $python_dir/bedgraph_rescale.py \
    -P $sample_dir/"$library_type"_minus_mask.bedgraph \
    -N $bg_minus \
    -A $annotation_gff \
    -S - \
    --align \
    --position $library_type \
    > $sample_dir/"$library_type"_scale_minus.txt

plus_scales=( $(tail -n 1 $sample_dir/"$library_type"_scale_plus.txt) )
minus_scales=( $(tail -n 1 $sample_dir/"$library_type"_scale_minus.txt) )
scale=$(printf "%.2f" $(echo "(${plus_scales[0]} + ${minus_scales[0]}) / 2" | bc -l))
bandwidth=$(printf "%.0f" $(echo "(${plus_scales[1]} + ${minus_scales[1]}) / 2" | bc -l))

if [ $(echo "$scale > $SCALE_CAP" | bc -l) -eq 1 ]
then    
    scale=$SCALE_CAP
    echo "hit scale cap of $SCALE_CAP"
fi

if [ $(echo "$bandwidth > $BANDWIDTH_CAP" | bc -l) -eq 1 ]
then
    bandwidth=$BANDWIDTH_CAP
    echo "hit bandwidth cap of $BANDWIDTH_CAP"
fi

python $python_dir/bedgraph_combine.py \
    -i $sample_dir/"$library_type"_plus_mask.bedgraph $bg_plus  \
    -s $scale -1 \
    -o $sample_dir/"$library_type"_plus_subtract.bedgraph

python $python_dir/bedgraph_combine.py \
    -i $sample_dir/"$library_type"_minus_mask.bedgraph $bg_minus \
    -s $scale -1 \
    -o $sample_dir/"$library_type"_minus_subtract.bedgraph

echo $scale > $sample_dir/final_scaling_factors.tab
rm $sample_dir/tmp.bedgraph

### PHASE 3.4: CONTINUOUS KERNEL DENSITY DISTRIBUTION ###
echo "#########################################################"
echo "### PHASE 3.4: CONTINUOUS KERNEL DENSITY DISTRIBUTION ###"
echo "#########################################################"
echo " "
# Takes the four subtractive BEDGRAPH files generated in PHASE 3 
# and smooths them using a continous kernel function.
# See the Python util "bedgraph_kernel_density.py" for full details.
# Defaults to a summed Laplace distribution smoothing

for strand in plus minus
do
    kernel_density_command="python \
        $python_dir/bedgraph_kernel_density.py \
        -B=$sample_dir/"$library_type"_"$strand"_subtract.bedgraph \
        -O=$sample_dir/"$library_type"_"$strand"_smooth.bedgraph \
        -L=$resource_dir/length.table \
        -K=$KERNEL \
        -H=$bandwidth \
        -S=3 \
        -D=3 \
        -P \
        -c $CPUS"
    echo $kernel_density_command
    eval $kernel_density_command
    if [ ! -f $sample_dir/"$library_type"_"$strand"_smooth.bedgraph ]
    then
        echo "ERROR: Failed to generate "$library_type"_"$strand"_smooth.bedgraph"
        exit 1
    fi
    feature_threshold_command="python \
    $python_dir/bedgraph_thresh_to_bed.py \
    -B $sample_dir/"$library_type"_"$strand"_smooth.bedgraph \
    -O $sample_dir/"$library_type"_"$strand"_features.bed \
    -L $resource_dir/length.table \
    -T 0 \
    -M 10 \
    -V sum \
    -S $strand"
    echo $feature_threshold_command
    eval $feature_threshold_command
    if [ ! -f $sample_dir/"$library_type"_"$strand"_features.bed ]
    then
        echo "ERROR: Failed to generate "$library_type"_"$strand"_features.bed"
        exit 1
    fi
done

# Merges all end features identified in PHASE 3.4 to a single BED file.
rm -f $sample_dir/end_features_temp.bed
touch $sample_dir/end_features_temp.bed
sed 's/thresh./'$library_type'.plus./' $sample_dir/"$library_type"_plus_features.bed \
    >> $sample_dir/end_features_temp.bed
sed 's/thresh./'$library_type'.minus./' $sample_dir/"$library_type"_minus_features.bed \
    >> $sample_dir/end_features_temp.bed

bedtools sort -i $sample_dir/end_features_temp.bed > $sample_dir/end_features.bed

python $python_dir/bed_find_peaks.py \
    -I $sample_dir/end_features.bed \
    -L $resource_dir/length.table \
    -O $sample_dir/$sample_name.rpm.features.bed \
    -BP $sample_dir/"$library_type"_plus_mask.bedgraph \
    -BM $sample_dir/"$library_type"_minus_mask.bedgraph \
    -V rpm

awk -F'[\t]' '{if ($5 >= .5){ print }}' $sample_dir/$sample_name.rpm.features.bed > $sample_dir/$sample_name.end_features.bed

rm $sample_dir/end_features_temp.bed \
    $sample_dir/end_features.bed \
    $sample_dir/$sample_name.rpm.features.bed \
    $sample_dir/"$library_type"_plus_features.bed \
    $sample_dir/"$library_type"_minus_features.bed
    $sample_dir/"$library_type"_plus_subtract.bedgraph \
    $sample_dir/"$library_type"_minus_subtract.bedgraph

for strand in plus minus
do
    rm $sample_dir/"$library_type"_"$strand".bedgraph
done

echo "Phase 3.4 complete."

echo "EndGraph complete!"


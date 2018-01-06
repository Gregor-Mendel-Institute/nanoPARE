#!/bin/bash

# EndGraph
# Identification of end features through subtractive kernel density estimation
# 1) Determines a scaling factor to adjust read depths of 5P/3P and BODY libraries
# 2) Signal smoothing by fitting a kernel to END - BODY read values
# 3) Convert continuous regions of positive signal to features in a BED file

################
# CONFIG SETUP #
################
# Storing all default global environment variables

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash_dir=$root_dir/scripts/bash_scripts
python_dir=$root_dir/scripts/python_scripts
r_dir=$root_dir/scripts/r_scripts
resource_dir=$root_dir/resources
temp_dir=$root_dir/temp
log_dir=$root_dir/log
results_dir=$root_dir/results

genome_fasta=$resource_dir/genome.fasta
annotation_gff=$resource_dir/annotation.gff
TSS_PLUS=$resource_dir/TSS_plus.bedgraph
TSS_MINUS=$resource_dir/TSS_minus.bedgraph
TES_PLUS=$resource_dir/TES_plus.bedgraph
TES_MINUS=$resource_dir/TES_minus.bedgraph
BODY_PLUS=$resource_dir/BODY_plus.bedgraph
BODY_MINUS=$resource_dir/BODY_minus.bedgraph
UUG_PLUS=
UUG_MINUS=

SAMPLE_NAME="sample"
LMOD=0
CPUS=1
SETUP=false
KERNEL='laplace'

############################
# READING THE COMMAND LINE #
############################
# Taking the default variables above, modifying them with the commandline 
# arguments (see read_cmdline.sh), and writing a config file

. $bash_dir/read_cmdline.sh
temp_dir=$temp_dir/$SAMPLE_NAME
mkdir -p $temp_dir

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
        $resource_dir
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
# where TSS or TES reads are likely methodological artifacts:
#   TSO oligos perform strand invasion at internal sites complementary
#     to their 3' ends (Tang et al. 2012, Nucleic Acids Research)
#   oligo(dT) priming results in internal priming artifacts at templated 
#     A-rich or purine-rich regions of transcripts
# Customize what sequences are masked by editing the file 
# /resources/mask_sequences.table

if [ -f $TSS_PLUS ] && [ -f $TSS_MINUS ]
then
    TSS="true"
    cp $TSS_PLUS $temp_dir/TSS_plus.bedgraph
    cp $TSS_MINUS $temp_dir/TSS_minus.bedgraph
else
    TSS=false
fi
if [ -f $TES_PLUS ] && [ -f $TES_MINUS ]
then
    TES="true"
    cp $TES_PLUS $temp_dir/TES_plus.bedgraph
    cp $TES_MINUS $temp_dir/TES_minus.bedgraph
else
    TES=false
fi
cp $BODY_PLUS $temp_dir/BODY_plus.bedgraph
if [ -f $BODY_MINUS ]
then
    cp $BODY_MINUS $temp_dir/BODY_minus.bedgraph
else
    cp $BODY_PLUS $temp_dir/BODY_minus.bedgraph
fi

if [ $TSS == "true" ]
then
    python $python_dir/bedgraph_mask.py \
        -P $temp_dir/TSS_plus.bedgraph \
        -M $temp_dir/TSS_minus.bedgraph \
        -PO $temp_dir/TSS_plus_mask.bedgraph \
        -MO $temp_dir/TSS_minus_mask.bedgraph \
        -U $resource_dir/TSS_mask_up.bed \
        -L $resource_dir/length.table
fi
if [ $TES == "true" ]
then
    python $python_dir/bedgraph_mask.py \
        -P $temp_dir/TES_plus.bedgraph \
        -M $temp_dir/TES_minus.bedgraph \
        -PO $temp_dir/TES_plus_mask.bedgraph \
        -MO $temp_dir/TES_minus_mask.bedgraph \
        -D $resource_dir/TES_mask_down.bed \
        -L $resource_dir/length.table
fi

if [ -z $UUG_PLUS ]
then
    echo "uuG files not provided. Skipping uuG masking."
else
    echo "uuG files: $UUG_PLUS $UUG_MINUS"
    python $python_dir/bedgraph_mask.py \
        -P $UUG_PLUS \
        -M $UUG_MINUS \
        -PO $temp_dir/uuG_plus_mask.bedgraph \
        -MO $temp_dir/uuG_minus_mask.bedgraph \
        -U $resource_dir/TSS_mask_up.bed \
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
#     Scaling factor to adjust for capture efficiency of TSS/TES/BODY reads
#
# Generates four scaled subtractive bedgraph files:
#     TSSplus  - BODYminus
#     TSSminus - BODYplus
#     TESplus  - BODYplus
#     TESminus - BODYminus

TSS_scale=1
TES_scale=1
SCALE_CAP=10
BANDWIDTH_CAP=30

if [ $TSS == "true" ]
then
    python $python_dir/bedgraph_rescale.py \
        -P $temp_dir/TSS_plus_mask.bedgraph \
        -N $temp_dir/BODY_minus.bedgraph \
        -A $annotation_gff \
        -S + \
        --position TSS \
        --align \
        > $temp_dir/TSS_scale_plus.txt
    
    python $python_dir/bedgraph_rescale.py \
        -P $temp_dir/TSS_minus_mask.bedgraph \
        -N $temp_dir/BODY_plus.bedgraph \
        -A $annotation_gff \
        -S - \
        --align \
        --position TSS \
        > $temp_dir/TSS_scale_minus.txt
    
    TSS_plus_scales=( $(tail -n 1 $temp_dir/TSS_scale_plus.txt) )
    TSS_minus_scales=( $(tail -n 1 $temp_dir/TSS_scale_minus.txt) )
    TSS_scale=$(printf "%.2f" $(echo "(${TSS_plus_scales[0]} + ${TSS_minus_scales[0]}) / 2" | bc -l))
    TSS_bandwidth=$(printf "%.0f" $(echo "(${TSS_plus_scales[1]} + ${TSS_minus_scales[1]}) / 2" | bc -l))
    
    if [ $(echo "$TSS_scale > $SCALE_CAP" | bc -l) -eq 1 ]
    then    
        TSS_scale=$SCALE_CAP
        echo "TSS hit scale cap of $SCALE_CAP"
    fi

    if [ $(echo "$TSS_BANDWIDTH > $BANDWIDTH_CAP" | bc -l) -eq 1 ]
    then
        TSS_bandwidth=$BANDWIDTH_CAP
        echo "TSS hit bandwidth cap of $BANDWIDTH_CAP"
    fi

    
    unionBedGraphs -i $temp_dir/TSS_plus_mask.bedgraph \
        $temp_dir/BODY_minus.bedgraph \
        > $temp_dir/tmp.bedgraph
    
    awk -v mult=$TSS_scale \
        '{printf $1"\t"$2"\t"$3"\t"$4*mult-$5"\n"}' $temp_dir/tmp.bedgraph \
        > $temp_dir/TSS_plus_subtract.bedgraph
    
    unionBedGraphs -i $temp_dir/TSS_minus_mask.bedgraph \
        $temp_dir/BODY_plus.bedgraph \
        > $temp_dir/tmp.bedgraph
    
    awk -v mult=$TSS_scale \
        '{printf $1"\t"$2"\t"$3"\t"$4*mult-$5"\n"}' $temp_dir/tmp.bedgraph \
        > $temp_dir/TSS_minus_subtract.bedgraph
fi

if [ $TES == "true" ]
then
    python $python_dir/bedgraph_rescale.py \
        -P $temp_dir/TES_plus_mask.bedgraph \
        -N $temp_dir/BODY_plus.bedgraph \
        -A $annotation_gff \
        -S + \
        --position TES \
        --align \
        > $temp_dir/TES_scale_plus.txt
        
    python $python_dir/bedgraph_rescale.py \
        -P $temp_dir/TES_minus_mask.bedgraph \
        -N $temp_dir/BODY_minus.bedgraph \
        -A $annotation_gff \
        -S - \
        --position TES \
        --align \
        > $temp_dir/TES_scale_minus.txt
    
    TES_plus_scales=( $(tail -n 1 $temp_dir/TES_scale_plus.txt) )
    TES_minus_scales=( $(tail -n 1 $temp_dir/TES_scale_minus.txt) )
    TES_scale=$(printf "%.2f" $(echo "(${TES_plus_scales[0]} + ${TES_minus_scales[0]}) / 2" | bc -l))
    TES_bandwidth=$(printf "%.0f" $(echo "(${TES_plus_scales[2]} + ${TES_minus_scales[2]}) / 2" | bc -l))
    
    if [ $(echo "$TES_scale > $SCALE_CAP" | bc -l) -eq 1 ]
    then
        TES_scale=$SCALE_CAP
        echo "TES hit scale cap of $SCALE_CAP"
    fi
    
    if [ $(echo "$TES_bandwidth > $BANDWIDTH_CAP" | bc -l) -eq 1 ]
    then
        TES_bandwidth=$BANDWIDTH_CAP
        echo "TES hit bandwidth cap of $BANDWIDTH_CAP"
    fi

    unionBedGraphs -i $temp_dir/TES_plus_mask.bedgraph \
        $temp_dir/BODY_plus.bedgraph \
        > $temp_dir/tmp.bedgraph

    awk -v mult=$TES_scale \
        '{printf $1"\t"$2"\t"$3"\t"$4*mult-$5"\n"}' $temp_dir/tmp.bedgraph \
        > $temp_dir/TES_plus_subtract.bedgraph

    unionBedGraphs -i $temp_dir/TES_minus_mask.bedgraph \
        $temp_dir/BODY_minus.bedgraph \
        > $temp_dir/tmp.bedgraph

    awk -v mult=$TES_scale \
        '{printf $1"\t"$2"\t"$3"\t"$4*mult-$5"\n"}' $temp_dir/tmp.bedgraph \
        > $temp_dir/TES_minus_subtract.bedgraph
fi

echo "$TSS_scale $TES_scale" > $temp_dir/final_scaling_factors.tab
rm $temp_dir/tmp.bedgraph

### PHASE 3.4: CONTINUOUS KERNEL DENSITY DISTRIBUTION ###
echo "#########################################################"
echo "### PHASE 3.4: CONTINUOUS KERNEL DENSITY DISTRIBUTION ###"
echo "#########################################################"
echo " "
# Takes the four subtractive BEDGRAPH files generated in PHASE 3 
# and smooths them using a continous kernel function.
# See the Python util "bedgraph_kernel_density.py" for full details.
# Defaults to a summed Laplace distribution smoothing

if [ $TSS == "true" ] && [ $TES == "true" ]
then
    readtypes=( TSS TES )
elif [ $TSS == "true" ]
then
    readtypes=( TSS )
else
    readtypes=( TES )
fi

for readtype in ${readtypes[@]}
do
    if [ $readtype == 'TSS' ]
    then
        bandwidth=$TSS_bandwidth
    else
        bandwidth=$TES_bandwidth
    fi
    for strand in plus minus
    do
        kernel_density_command="python \
            $python_dir/bedgraph_kernel_density.py \
            -B=$temp_dir/"$readtype"_"$strand"_subtract.bedgraph \
            -O=$temp_dir/"$readtype"_"$strand"_smooth.bedgraph \
            -L=$resource_dir/length.table \
            -K=$KERNEL \
            -H=$bandwidth \
            -S=3 \
            -D=3 \
            -P \
            -c $CPUS"
        echo $kernel_density_command
        eval $kernel_density_command
        if [ ! -f $temp_dir/"$readtype"_"$strand"_smooth.bedgraph ]
        then
            echo "ERROR: Failed to generate "$readtype"_"$strand"_smooth.bedgraph"
            exit 1
        fi
        feature_threshold_command="python \
        $python_dir/bedgraph_thresh_to_bed.py \
        -B $temp_dir/"$readtype"_"$strand"_smooth.bedgraph \
        -O $temp_dir/"$readtype"_"$strand"_features.bed \
        -L $resource_dir/length.table \
        -T 0 \
        -M 10 \
        -V sum \
        -S $strand"
        echo $feature_threshold_command
        eval $feature_threshold_command
        if [ ! -f $temp_dir/"$readtype"_"$strand"_features.bed ]
        then
            echo "ERROR: Failed to generate "$readtype"_"$strand"_features.bed"
            exit 1
        fi
    done
done

# Merges all end features identified in PHASE 3.4 to a single BED file.
touch $temp_dir/end_features_temp.bed
if [ $TSS == "true" ]
then
    sed 's/thresh./TSS.plus./' $temp_dir/TSS_plus_features.bed \
        >> $temp_dir/end_features_temp.bed
    sed 's/thresh./TSS.minus./' $temp_dir/TSS_minus_features.bed \
        >> $temp_dir/end_features_temp.bed
fi
if [ $TES == "true" ]
then
    sed 's/thresh./TES.plus./' $temp_dir/TES_plus_features.bed \
        >> $temp_dir/end_features_temp.bed
    sed 's/thresh./TES.minus./' $temp_dir/TES_minus_features.bed \
        >> $temp_dir/end_features_temp.bed
fi
bedtools sort -i $temp_dir/end_features_temp.bed > $temp_dir/end_features.bed

python $python_dir/bed_find_peaks.py \
    -I $temp_dir/end_features.bed \
    -L $resource_dir/length.table \
    -O $temp_dir/$SAMPLE_NAME.rpm.features.bed \
    -BP $temp_dir/TSS_plus_mask.bedgraph \
    -BM $temp_dir/TSS_minus_mask.bedgraph \
    -V rpm

awk -F'[\t]' '{if ($5 >= .5){ print }}' $temp_dir/$SAMPLE_NAME.rpm.features.bed > $temp_dir/$SAMPLE_NAME.end_features.bed

rm $temp_dir/end_features_temp.bed \
    $temp_dir/end_features.bed \
    $temp_dir/$SAMPLE_NAME.rpm.features.bed \
    $temp_dir/TSS_plus_features.bed \
    $temp_dir/TSS_minus_features.bed \
    $temp_dir/TES_plus_features.bed \
    $temp_dir/TES_minus_features.bed \

rm $temp_dir/TSS_plus_subtract.bedgraph \
    $temp_dir/TSS_minus_subtract.bedgraph \
    $temp_dir/TES_plus_subtract.bedgraph \
    $temp_dir/TSS_minus_subtract.bedgraph

for readtype in ${readtypes[@]}
do
    for strand in plus minus
    do
        rm $temp_dir/"$readtype"_"$strand".bedgraph
    done
done

echo "Phase 3.4 complete."


echo "PHASE 3 complete!"



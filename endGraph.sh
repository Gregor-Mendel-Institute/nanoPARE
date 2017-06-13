#!/bin/bash

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

SAMPLE_NAME="sample"
LMOD=0
SETUP=false
ITERATIONS=5
KERNEL='laplace'

############################
# READING THE COMMAND LINE #
############################
# Taking the default variables above, modifying them with the commandline arguments (see read_cmdline.sh), and writing a config file

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
REQUIRED_MODULES=( --bedtools --r )
. $bash_dir/load_modules.sh
echo " "

if [ $SETUP = true ]
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

    # A GFF file of features centered around TSS, and a separate GFF file for TES.
    TSS_GFF_command="python \
    $python_dir/gff_feature_ends.py \
    -I=$annotation_gff \
    -O=$resource_dir/TSS_sites.gff \
    -E=TSS \
    -D=500 \
    -F=gene"
    echo "$TSS_GFF_command"
    eval "$TSS_GFF_command"

    TES_GFF_command="python \
    $python_dir/gff_feature_ends.py \
    -I=$annotation_gff \
    -O=$resource_dir/TES_sites.gff \
    -E=TES \
    -D=500 \
    -F=gene"
    echo "$TES_GFF_command"
    eval "$TES_GFF_command"
    echo "TSS and TES GFF files generated."

    python $python_dir/fasta_sequence_search.py $genome_fasta $resource_dir/mask_sequences.table $resource_dir
    echo "Masking BED files generated."
    echo "Setup complete."
    exit 0
fi


### STEP 1: BEDGRAPH ARTIFACT MASKING ###
echo "#########################################"
echo "### STEP 1: BEDGRAPH ARTIFACT MASKING ###"
echo "#########################################"
echo " "
# Generates a list of BED features indicating sites in the genome where TSS or TES reads are likely methodological artifacts:
#     TSO oligos perform strand invasion at internal sites complementary to their 3' ends (Tang et al. 2012, Nucleic Acids Research)
#     oligo(dT) priming results in internal priming artifacts at templated A-rich regions of transcripts
# Customize what sequences are masked by editing /resources/mask_sequences.table

if [ -f $TSS_PLUS ] && [ -f $TSS_MINUS ]
then
    TSS=true
    cp $TSS_PLUS $temp_dir/TSS_plus.bedgraph
    cp $TSS_MINUS $temp_dir/TSS_minus.bedgraph
else
    TSS=false
fi
if [ -f $TES_PLUS ] && [ -f $TES_MINUS ]
then
    TES=true
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

for strand in plus minus
do
    if [ $TSS = true ]
    then
        python $python_dir/bedgraph_mask.py -B=$temp_dir/TSS_"$strand".bedgraph -O=$temp_dir/TSS_"$strand"_mask.bedgraph -L=$resource_dir/length.table -S=$strand -U=$resource_dir/TSS_mask_up.bed 
    fi
    if [ $TES = true ]
    then
        python $python_dir/bedgraph_mask.py -B=$temp_dir/TES_"$strand".bedgraph -O=$temp_dir/TES_"$strand"_mask.bedgraph -L=$resource_dir/length.table -S=$strand -D=$resource_dir/TES_mask_down.bed
    fi
done

echo "Step 1 complete."

### STEP 2: PARAMETER ESTIMATION WITH SUBTRACTIVE BEDGRAPH FILES ###
echo "####################################################################"
echo "### STEP 2: PARAMETER ESTIMATION WITH SUBTRACTIVE BEDGRAPH FILES ###"
echo "####################################################################"
echo " "
# Generates four subtractive bedgraph files based on theoretical read distribution:
#     TSSplus  - BODYminus
#     TSSminus - BODYplus
#     TESplus  - BODYplus
#     TESminus - BODYminus
# 
# Performs a meta-analysis on annotated TSS/TES sites to determine 2 parameters:
#     Optimal bandwidth for kernel density estimation
#     Scaling factor to adjust for varying capture efficiency of TSS, TES, and BODY reads
# 
# Rebuilds subtractive bedgraph files with the appropriate scaling factor

TSS_scale=1
TES_scale=1
TSS_bandwidth=50
TES_bandwidth=50

if [ $TSS = true ]
then
    unionBedGraphs -i $temp_dir/TSS_plus_mask.bedgraph $temp_dir/BODY_minus.bedgraph > $temp_dir/tmp.bedgraph
    awk -v mult=$TSS_scale '{printf $1"\t"$2"\t"$3"\t"$4/mult-$5*mult"\n"}' $temp_dir/tmp.bedgraph > $temp_dir/TSS_plus_subtract.bedgraph
    unionBedGraphs -i $temp_dir/TSS_minus_mask.bedgraph $temp_dir/BODY_plus.bedgraph > $temp_dir/tmp.bedgraph
    awk -v mult=$TSS_scale '{printf $1"\t"$2"\t"$3"\t"$4/mult-$5*mult"\n"}' $temp_dir/tmp.bedgraph > $temp_dir/TSS_minus_subtract.bedgraph
fi
if [ $TES = true ]
then
    unionBedGraphs -i $temp_dir/TES_plus_mask.bedgraph $temp_dir/BODY_plus.bedgraph > $temp_dir/tmp.bedgraph
    awk -v mult=$TES_scale '{printf $1"\t"$2"\t"$3"\t"$4/mult-$5*mult"\n"}' $temp_dir/tmp.bedgraph > $temp_dir/TES_plus_subtract.bedgraph
    unionBedGraphs -i $temp_dir/TES_minus_mask.bedgraph $temp_dir/BODY_minus.bedgraph > $temp_dir/tmp.bedgraph
    awk -v mult=$TES_scale '{printf $1"\t"$2"\t"$3"\t"$4/mult-$5*mult"\n"}' $temp_dir/tmp.bedgraph > $temp_dir/TES_minus_subtract.bedgraph
fi

if [ $ITERATIONS -ge 1 ]
then
    echo "Read scaling optimization:"
    for ((run=1;run<=$ITERATIONS;run++))
    do
        if [ $TSS = true ]
        then
            python $python_dir/bedgraph_metaplot.py -G=$resource_dir/TSS_sites.gff -B=$temp_dir/TSS_plus_subtract.bedgraph -S=plus -L=$resource_dir/length.table > $temp_dir/TSS_start_plus.meta
            python $python_dir/bedgraph_metaplot.py -G=$resource_dir/TSS_sites.gff -B=$temp_dir/TSS_minus_subtract.bedgraph -S=minus -L=$resource_dir/length.table > $temp_dir/TSS_start_minus.meta
            python $python_dir/bedgraph_metaplot.py -G=$resource_dir/TES_sites.gff -B=$temp_dir/TSS_plus_subtract.bedgraph -S=plus -L=$resource_dir/length.table > $temp_dir/TSS_end_plus.meta
            python $python_dir/bedgraph_metaplot.py -G=$resource_dir/TES_sites.gff -B=$temp_dir/TSS_minus_subtract.bedgraph -S=minus -L=$resource_dir/length.table > $temp_dir/TSS_end_minus.meta

            paste $temp_dir/TSS_start_plus.meta $temp_dir/TSS_start_minus.meta $temp_dir/TSS_end_plus.meta $temp_dir/TSS_end_minus.meta | awk '{printf "%.8f %.8f\n", ($1+$2)/2, ($3+$4)/2}' > $temp_dir/TSS_metaplot.tab
            rm $temp_dir/TSS_start_plus.meta $temp_dir/TSS_start_minus.meta $temp_dir/TSS_end_plus.meta $temp_dir/TSS_end_minus.meta

            cp $temp_dir/TSS_metaplot.tab $temp_dir/TSS_metaplot_$run.tab
            Rscript $r_dir/TSS_TES_scaling_factors.R $temp_dir/TSS_metaplot_$run.tab > $temp_dir/TSS_scaling_factors_$run.tab
            TSS_scales=$(echo -n $(for ((i=1;i<=$run;i++)); do echo $(cut -d ' ' -f 1 $temp_dir/TSS_scaling_factors_$i.tab); done | tr '\n' ' ') | tr ' ' '*')
            TSS_scale=$(echo $TSS_scales | bc)
            TSS_bandwidth=$(cut -d ' ' -f 2 $temp_dir/TSS_scaling_factors_$run.tab)
            echo "TSS_run_$run $TSS_scale $TSS_bandwidth"
            unionBedGraphs -i $temp_dir/TSS_plus_mask.bedgraph $temp_dir/BODY_minus.bedgraph > $temp_dir/tmp.bedgraph
            awk -v mult=$TSS_scale '{printf $1"\t"$2"\t"$3"\t"$4/(1+((mult-1)/2))-$5*(1+((mult-1)/2))"\n"}' $temp_dir/tmp.bedgraph > $temp_dir/TSS_plus_subtract.bedgraph
            unionBedGraphs -i $temp_dir/TSS_minus_mask.bedgraph $temp_dir/BODY_plus.bedgraph > $temp_dir/tmp.bedgraph
            awk -v mult=$TSS_scale '{printf $1"\t"$2"\t"$3"\t"$4/(1+((mult-1)/2))-$5*(1+((mult-1)/2))"\n"}' $temp_dir/tmp.bedgraph > $temp_dir/TSS_minus_subtract.bedgraph
        fi
        if [ $TES = true ]
        then
            python $python_dir/bedgraph_metaplot.py -G=$resource_dir/TSS_sites.gff -B=$temp_dir/TES_plus_subtract.bedgraph -S=plus -L=$resource_dir/length.table > $temp_dir/TES_start_plus.meta
            python $python_dir/bedgraph_metaplot.py -G=$resource_dir/TSS_sites.gff -B=$temp_dir/TES_minus_subtract.bedgraph -S=minus -L=$resource_dir/length.table > $temp_dir/TES_start_minus.meta
            python $python_dir/bedgraph_metaplot.py -G=$resource_dir/TES_sites.gff -B=$temp_dir/TES_plus_subtract.bedgraph -S=plus -L=$resource_dir/length.table > $temp_dir/TES_end_plus.meta
            python $python_dir/bedgraph_metaplot.py -G=$resource_dir/TES_sites.gff -B=$temp_dir/TES_minus_subtract.bedgraph -S=minus -L=$resource_dir/length.table > $temp_dir/TES_end_minus.meta

            paste $temp_dir/TES_start_plus.meta $temp_dir/TES_start_minus.meta $temp_dir/TES_end_plus.meta $temp_dir/TES_end_minus.meta | awk '{printf "%.8f %.8f\n", ($1+$2)/2, ($3+$4)/2}' > $temp_dir/TES_metaplot.tab
            rm $temp_dir/TES_start_plus.meta $temp_dir/TES_start_minus.meta $temp_dir/TES_end_plus.meta $temp_dir/TES_end_minus.meta

            cp $temp_dir/TES_metaplot.tab $temp_dir/TES_metaplot_$run.tab
            Rscript $r_dir/TSS_TES_scaling_factors.R $temp_dir/TSS_metaplot_$run.tab > $temp_dir/TES_scaling_factors_$run.tab
            TES_scales=$(echo -n $(for ((i=1;i<=$run;i++)); do echo $(cut -d ' ' -f 1 $temp_dir/TES_scaling_factors_$i.tab); done | tr '\n' ' ') | tr ' ' '*')
            TES_scale=$(echo $TES_scales | bc)
            TES_bandwidth=$(cut -d ' ' -f 2 $temp_dir/TES_scaling_factors_$run.tab)
            echo "TES_run_$run $TES_scale $TES_bandwidth"
            unionBedGraphs -i $temp_dir/TES_plus_mask.bedgraph $temp_dir/BODY_minus.bedgraph > $temp_dir/tmp.bedgraph
            awk -v mult=$TES_scale '{printf $1"\t"$2"\t"$3"\t"$4/(1+((mult-1)/2))-$5*(1+((mult-1)/2))"\n"}' $temp_dir/tmp.bedgraph > $temp_dir/TES_plus_subtract.bedgraph
            unionBedGraphs -i $temp_dir/TES_minus_mask.bedgraph $temp_dir/BODY_plus.bedgraph > $temp_dir/tmp.bedgraph
            awk -v mult=$TES_scale '{printf $1"\t"$2"\t"$3"\t"$4/(1+((mult-1)/2))-$5*(1+((mult-1)/2))"\n"}' $temp_dir/tmp.bedgraph > $temp_dir/TES_minus_subtract.bedgraph
        fi
    done
fi
rm $temp_dir/tmp.bedgraph

### PHASE 3.4: CONTINUOUS KERNEL DENSITY DISTRIBUTION ###
echo "#########################################################"
echo "### PHASE 3.4: CONTINUOUS KERNEL DENSITY DISTRIBUTION ###"
echo "#########################################################"
echo " "
# Takes the four subtractive BEDGRAPH files generated in PHASE 3.3 and smooths them using a continous kernel function.
# See the Python util "bedgraph_kernel_density.py" for full details.
# Defaults to a summed Laplace distribution smoothing

if [ $TSS = true ] && [ $TES = true ]
then
    readtypes=( TSS TES )
elif [ $TSS = true ]
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
        -S=5 \
        -D=3 \
        -P=True"
        echo $kernel_density_command
        eval $kernel_density_command
        if [ ! -f $temp_dir/"$readtype"_"$strand"_smooth.bedgraph ]
        then
            echo "ERROR: Failed to generate "$readtype"_"$strand"_smooth.bedgraph"
            exit 1
        fi
        feature_threshold_command="python \
        $python_dir/bedgraph_thresh_to_bed.py \
        -B=$temp_dir/"$readtype"_"$strand"_smooth.bedgraph \
        -O=$temp_dir/"$readtype"_"$strand"_features.bed \
        -L=$resource_dir/length.table \
        -T=0 \
        -M=10 \
        -V=sum \
        -S=$strand"
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
if [ $TSS = true ]
then
    sed 's/thresh./TSS.plus./' $temp_dir/TSS_plus_features.bed >> $temp_dir/end_features_temp.bed
    sed 's/thresh./TSS.minus./' $temp_dir/TSS_minus_features.bed >> $temp_dir/end_features_temp.bed
    Rscript $r_dir/TSS_TES_metaplot.R $temp_dir/TSS_metaplot.tab
fi
if [ $TES = true ]
then
    sed 's/thresh./TES.plus./' $temp_dir/TES_plus_features.bed >> $temp_dir/end_features_temp.bed
    sed 's/thresh./TES.minus./' $temp_dir/TES_minus_features.bed >> $temp_dir/end_features_temp.bed
    Rscript $r_dir/TSS_TES_metaplot.R $temp_dir/TES_metaplot.tab
fi
bedtools sort -i $temp_dir/end_features_temp.bed > $temp_dir/end_features.bed
rm $temp_dir/end_features_temp.bed $temp_dir/TSS_plus_features.bed $temp_dir/TSS_minus_features.bed $temp_dir/TES_plus_features.bed $temp_dir/TES_minus_features.bed
# rm $temp_dir/TSS_plus_subtract.bedgraph $temp_dir/TSS_minus_subtract.bedgraph $temp_dir/TES_plus_subtract.bedgraph $temp_dir/TSS_minus_subtract.bedgraph
for readtype in ${readtypes[@]}
do
    for strand in plus minus
    do
    rm $temp_dir/"$readtype"_"$strand".bedgraph
    done
done

echo "Phase 3.4 complete."
echo "PHASE 3 complete!"



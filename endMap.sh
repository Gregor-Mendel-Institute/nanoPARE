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
    reference_table=$resource_dir/reference_table_EndMap.txt
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


############################
# READING THE COMMAND LINE #
############################
# Taking the default variables above, modifying them with the commandline 
# arguments (see read_cmdline.sh), and writing a config file

. $bash_dir/read_cmdline.sh

echo "##############"
echo "### ENDMAP ###"
echo "##############"
echo " "
echo "Config settings:"
. $bash_dir/list_settings.sh


# Environment modules to load with Lmod (if option --lmod is passed)
REQUIRED_MODULES=( --bedtools --cutadapt --rna-star --samtools )
. $bash_dir/load_modules.sh
echo " "

OPTIONS_star_global=$(cat $resource_dir/OPTIONS_star_global)
input_mapper=$(sed -n "$JOB_NUMBER"p $reference_table) #read mapping file
input_array=($input_mapper)

line_number=${input_array[0]}  # Line number of reference table (should match job number)
fastq_dir=${input_array[1]}    # Directory containing the fastq file
input_fastq=${input_array[2]}  # FASTQ filename(s), comma-separated
SAMPLE_NAME=${input_array[3]}  # Name of the sample
library_type=${input_array[4]} # Options: BODY, TSS. Type of the FASTQ file
read_type=${input_array[5]}    # Options: SE, PE, for single end or paired end libraries
adapter_str=${input_array[6]}  # comma-separated list of sequences to trim from the 3' end of reads
RAM=30

mkdir -p $temp_dir
sample_dir=$temp_dir/$SAMPLE_NAME
rm -rf $sample_dir
mkdir -p $sample_dir

if [ $read_type == "SE" ]
then
    if [ $library_type == "TSS" ]
    then
        star_params_allreads="$OPTIONS_star_global \
        --alignEndsType Local \
        --outSJfilterOverhangMin -1 10 10 10 \
        --outSJfilterCountTotalMin -1 2 2 2 \
        --outSJfilterCountUniqueMin -1 2 2 2 \
        --twopassMode Basic \
        --outFilterMismatchNoverLmax 0.05 \
        --outFilterMatchNminOverLread 0.9 \
        --readFilesIn $sample_dir/TSO_single_cleaned.fastq \
        --outFileNamePrefix $sample_dir/TSSstar_single/"
    else
        star_params_allreads="$OPTIONS_star_global \
        --alignEndsType EndToEnd \
        --outSJfilterOverhangMin -1 10 10 10 \
        --outSJfilterCountTotalMin -1 2 2 2 \
        --outSJfilterCountUniqueMin -1 2 2 2 \
        --twopassMode Basic \
        --outFileNamePrefix $sample_dir/star/ \
        --readFilesIn $sample_dir/"$sample_name"_cleaned.1.fastq"
    fi
elif [ $read_type == "PE" ]
    star_params_allreads="$OPTIONS_star_global \
    --alignEndsType EndToEnd \
    --outSJfilterOverhangMin -1 10 10 10 \
    --outSJfilterCountTotalMin -1 2 2 2 \
    --outSJfilterCountUniqueMin -1 2 2 2 \
    --twopassMode Basic \
    --outFileNamePrefix $sample_dir/star/ \
    --readFilesIn $sample_dir/"$sample_name"_cleaned.1.fastq $sample_dir/"$sample_name"_cleaned.2.fastq"
fi

#########################

### EXECUTE COMMANDS ###
echo "### PHASE 1: STAGE-SPECIFIC FASTQ FILE ASSEMBLY AND ADAPTER TRIMMING ###"
# Parse $input_fastq to determine:
# 1) Whether 1 or 2 files were specified in the reference.table
# 2) Whether the file(s) is/are compressed with gzip
# 3) Whether to split the file (a single PE FASTQ)
if [[ $input_fastq = *,* ]]
then
    if [ $read_type == "SE" ]
    then 
        echo "ERROR: Multiple FASTQ files provided for single-end experiment"
        exit 1
    fi
    echo "Splitting input_fastq into multiple files"
    fastq_split=$(echo $input_fastq | tr ',' ' ')
    fastq_array=($fastq_split)
    fastq_file_number=${#fastq_array[@]}
    if [ $fastq_file_number -gt 2 ]
    then
        echo "ERROR: Too many FASTQ files provided."
        echo "For paired-end experiments, please provide either one integrated FASTQ file,"
        echo "or a comma-separated pair of FASTQ files, one for each mate."
        exit 1
    fi
else
    fastq_array=($input_fastq)
    fastq_file_number=${#fastq_array[@]}
fi

# Gunzip any gzipped FASTQ files and move the files to the temp directory
mate_number=1
for fastq_file in "${fastq_array[@]}"
do
    if [[ $fastq_file = *.gz ]]
    then
        echo "Unzipping fastq file $fastq_file to the sample temp directory"
        gunzip -k -c $fastq_dir/$fastq_file > $sample_dir/"$sample_name"."$mate_number".fastq
    else
        echo "Copying fastq file $fastq_file to the sample temp directory"
        cat $fastq_dir/$fastq_file > $sample_dir/"$sample_name"."$mate_number".fastq
    fi
    ((mate_number+=1))
done

# If paired-end reads are specified and a single FASTQ file was provided, split the FASTQ into mate pairs 1 and 2
if [[ $read_type == "PE"* ]] && [ $fastq_file_number -eq 1 ]
then
    echo "Splitting paired-end FASTQ file..."
    mv $sample_dir/"$sample_name".1.fastq $sample_dir/"$sample_name".merged.fastq
    python $python_dir/fastq_split_paired_end.py $sample_dir/"$sample_name".merged.fastq $sample_dir/"$sample_name".1.fastq $sample_dir/"$sample_name".2.fastq
    rm $sample_dir/"$sample_name".merged.fastq
fi

pair_one=$sample_dir/$sample_name.1.fastq
pair_two=$sample_dir/$sample_name.2.fastq
touch $pair_one
touch $pair_two
read1_lines=$(wc -l < $pair_one)
read2_lines=$(wc -l < $pair_two)

if [ $read_type == "SE" ] && [ $read2_lines -gt 0 ]
then
    echo "ERROR: Non-empty mate pair 2 in SE experiment."
    exit 1
fi
if [[ $read_type == "PE"* ]] && [ ! $read1_lines -eq $read1_lines ]
then
    echo "ERROR: paired-end read numbers unequal."
    exit 1
fi

echo " "
echo "############################"
echo "Finished parsing FASTQ files"
echo "Experiment type: $read_type"
echo "Mate 1 reads: $(($read1_lines / 4))"
echo "Mate 2 reads: $(($read2_lines / 4))"
echo "############################"
echo " "
echo "Trimming transposase adapters with cutadapt:"

cd $sample_dir

# Trim 0, 1, or 2 provided adapter sequences

keep_untrimmed="true"

if [[ $library_type == "TSS" ]]
then
    adapters=( $(echo $adapter_str | tr "," " ") )
    number_of_adapters=${#adapters[@]}
    echo "Adapters: ${adapters[@]}"
    if [[ $number_of_adapters -eq 2 ]]
    then
        first_trim_command="cutadapt \
            -a ${adapters[0]} \
            -o "$SAMPLE_NAME"_trim1.fastq \
            --untrimmed-output "$SAMPLE_NAME"_untrimmed1.fastq \
            "$SAMPLE_NAME".fastq"

        second_trim_command="cutadapt \
            -a ${adapters[1]} \
            -o "$SAMPLE_NAME"_trim2.fastq \
            "$SAMPLE_NAME"_untrimmed1.fastq"

        echo $first_trim_command
        eval $first_trim_command
        echo $second_trim_command
        eval $second_trim_command

        cat "$SAMPLE_NAME"_trim1.fastq "$SAMPLE_NAME"_trim2.fastq > "$SAMPLE_NAME"_adaptertrim.fastq
        rm -f "$SAMPLE_NAME"_trim1.fastq "$SAMPLE_NAME"_trim2.fastq "$SAMPLE_NAME"_untrimmed1.fastq

    else
        if [[ $number_of_adapters -eq 1 ]]
        then
            trim_command="cutadapt \
                -a ${adapters[0]} \
                -o "$SAMPLE_NAME"_adaptertrim.fastq \
                --untrimmed-output "$SAMPLE_NAME"_untrimmed1.fastq
                "$SAMPLE_NAME".fastq"

            echo $trim_command
            eval $trim_command

        else
            if [[ $keep_untrimmed == "true" ]]
            then
                cat "$SAMPLE_NAME".fastq "$SAMPLE_NAME"_untrimmed1.fastq > "$SAMPLE_NAME"_adaptertrim.fastq
            else
                cat "$SAMPLE_NAME".fastq > "$SAMPLE_NAME"_adaptertrim.fastq
            fi
        fi
    fi
    cat "$SAMPLE_NAME"_adaptertrim.fastq > "$SAMPLE_NAME"_trimmed.fastq
    python $python_dir/fastq_drop_short_reads.py . "$SAMPLE_NAME"_trimmed.fastq "$SAMPLE_NAME"_cleaned.fastq 16

    echo "Removing trim intermediates"
    rm -f "$SAMPLE_NAME"_trim1.fastq "$SAMPLE_NAME"_adaptertrim.fastq "$SAMPLE_NAME"_trimmed.fastq
    echo "Adapter trimming complete."

    echo "Filtering out low-complexity reads..."
    python $python_dir/fastq_complexity_filter.py $sample_dir "$SAMPLE_NAME"_cleaned.fastq "$SAMPLE_NAME"_filtered.fastq 0.15
    rm -f "$SAMPLE_NAME"_cleaned.fastq
else
    #TODO: Improve adapter trimming behavior for BODY reads
    TN5_1='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
    TN5_2='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
    TN5_1rc=$(python $python_dir/stdin_rc.py $TN5_1 RC)
    TN5_2rc=$(python $python_dir/stdin_rc.py $TN5_2 RC)
    echo "tn5.1: $TN5_1 (reverse complement $TN5_1rc)"
    echo "tn5.2: $TN5_2 (reverse complement $TN5_2rc)"

    echo "cutadapt -a $TN5_2rc -o "$sample_name"_adaptertrim.1.fastq "$sample_name".1.fastq #rc(Tn5.2)"
    eval "cutadapt -a $TN5_2rc -o "$sample_name"_adaptertrim.1.fastq "$sample_name".1.fastq #rc(Tn5.2)"

    echo "cutadapt -a $TN5_1rc -o "$sample_name"_adaptertrim.2.fastq "$sample_name".2.fastq #rc(Tn5.1)"
    eval "cutadapt -a $TN5_1rc -o "$sample_name"_adaptertrim.2.fastq "$sample_name".2.fastq #rc(Tn5.1)"

    ###########
    if [[ $read_type == "PE"* ]]
    then
        echo "Cleaning reads with fastq_drop_short_pairs.py"
        echo "python $python_dir/fastq_drop_short_pairs.py $sample_dir "$sample_name"_adaptertrim.1.fastq "$sample_name"_adaptertrim.2.fastq "$sample_name"_cleaned.1.fastq "$sample_name"_cleaned.2.fastq $minimum_pairlength"
        eval "python $python_dir/fastq_drop_short_pairs.py $sample_dir "$sample_name"_adaptertrim.1.fastq "$sample_name"_adaptertrim.2.fastq "$sample_name"_cleaned.1.fastq "$sample_name"_cleaned.2.fastq $minimum_pairlength"
        rm "$sample_name"_adaptertrim.1.fastq "$sample_name"_adaptertrim.2.fastq "$sample_name".1.fastq "$sample_name".2.fastq
    else
        echo "Cleaning reads with fastq_drop_short_reads.py"
        echo "python $python_dir/fastq_drop_short_reads.py $sample_dir "$sample_name"_adaptertrim.1.fastq "$sample_name"_cleaned.1.fastq $minimum_readlength"
        eval "python $python_dir/fastq_drop_short_reads.py $sample_dir "$sample_name"_adaptertrim.1.fastq "$sample_name"_cleaned.1.fastq $minimum_readlength"
        rm "$sample_name"_adaptertrim.1.fastq "$sample_name"_adaptertrim.2.fastq "$sample_name".1.fastq "$sample_name".2.fastq
    fi
    echo "Adapter trimming complete."
fi

echo "### PHASE 2: MAPPING READS TO THE GENOME WITH STAR ###"
mkdir -p $sample_dir/star
echo "STAR $star_params_allreads >& $SAMPLE_NAME.star.log"
eval "STAR $star_params_allreads >& $SAMPLE_NAME.star.log"

rm "$SAMPLE_NAME".fastq

echo Generating bedgraph files...
samtools view -h star/Aligned.out.bam > "$SAMPLE_NAME".sam

if [[ $library_type == "TSS" ]]
then
    python $python_dir/sam_calculate_coverages.py \
        -S "$SAMPLE_NAME".5p \
        -I "$SAMPLE_NAME".sam \
        -R TSS \
        -F $genome_fasta \
        --softclip_type 5p \
        --untemp_out A C G T \
        --secondary \
        --allow_naive
else
    python $python_dir/sam_calculate_coverages.py \
        -S "$SAMPLE_NAME".5p \
        -I "$SAMPLE_NAME".sam \
        -R $library_type \
        -F $genome_fasta \
        --secondary \
        --allow_naive
fi

rm "$SAMPLE_NAME".sam 

bedtools sort -i "$SAMPLE_NAME".5p_"$library_type"_plus.bedgraph > "$SAMPLE_NAME"_plus.5p.bedgraph
bedtools sort -i "$SAMPLE_NAME".5p_"$library_type"_minus.bedgraph > "$SAMPLE_NAME"_minus.5p.bedgraph

if [[ $library_type == "TSS" ]]
then
    bedtools sort -i "$SAMPLE_NAME".5p_"$library_type"_plus_untemp.bedgraph | \
        awk '{ printf $1"\t"$2"\t"$3"\t"$6"\n }' > "$SAMPLE_NAME"_plus.uG.bedgraph
    bedtools sort -i "$SAMPLE_NAME".5p_"$library_type"_minus_untemp.bedgraph | \
        awk '{ printf $1"\t"$2"\t"$3"\t"$6"\n }' > "$SAMPLE_NAME"_minus.uG.bedgraph
    # Calculate transcript_level coverage by parsing the genome coverage values
    python $python_dir/bedgraph_genome_to_transcripts.py \
        --subset $annotation_subset \
        --output "$SAMPLE_NAME"_transcript.bedgraph \
        "$SAMPLE_NAME"_plus.5p.bedgraph \
        "$SAMPLE_NAME"_minus.5p.bedgraph \
        $annotation_gff \
        $genome_fasta
else
    # Calculate transcript_level coverage by parsing the genome coverage values
    python $python_dir/bedgraph_genome_to_transcripts.py \
        --subset $annotation_subset \
        --output "$SAMPLE_NAME"_transcript.bedgraph \
        "$SAMPLE_NAME"_minus.5p.bedgraph \
        "$SAMPLE_NAME"_plus.5p.bedgraph \
        $annotation_gff \
        $genome_fasta
    
fi

rm "$SAMPLE_NAME".5p_"$library_type"_plus.bedgraph \
    "$SAMPLE_NAME".5p_"$library_type"_minus.bedgraph

# rm "$SAMPLE_NAME"*coverage.bedgraph

echo "Moving final files to results folder..."
mkdir -p $results_dir/EndMap/$SAMPLE_NAME
cp $sample_dir/*.bedgraph $results_dir/EndMap/$SAMPLE_NAME/
cat $sample_dir/comptable.tsv > $results_dir/EndMap/$SAMPLE_NAME/"$SAMPLE_NAME".comptable.tsv

echo Pipeline complete!

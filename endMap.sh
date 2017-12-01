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

input_mapper=$(sed -n "$JOB_NUMBER"p $reference_table) #read mapping file
input_array=($input_mapper)

line_number=${input_array[0]}  # Line number of reference table (should match job number)
sample_dir=${input_array[1]}   # Directory containing the fastq file
adapter_str=${input_array[2]}     # comma-separated list of sequences to trim from the 3' end of reads
SAMPLE_NAME=$(basename $sample_dir)

mkdir -p $temp_dir
temp_dir_s=$temp_dir/$SAMPLE_NAME
rm -rf $temp_dir_s

star_params_allreads="--runMode alignReads \
--runThreadN $CPUS \
--runRNGseed 12345 \
--genomeDir $genome_index_dir \
--readFilesIn $sample_dir/"$SAMPLE_NAME"_filtered.fastq \
--outFileNamePrefix $sample_dir/star/ \
--limitBAMsortRAM 30000000000 \
--alignEndsType Local \
--outTmpDir $temp_dir_s \
--alignIntronMax 5000 \
--alignSJDBoverhangMin 1 \
--outReadsUnmapped Fastx \
--outSAMtype BAM Unsorted \
--outSAMprimaryFlag AllBestScore \
--outSAMmultNmax 100 \
--outSAMattributes NH HI AS nM NM MD jM jI XS \
--outFilterMultimapNmax 100 \
--outFilterMismatchNmax 2 \
--outFilterMatchNmin 16 \
--outFilterMatchNminOverLread 0.9 \
--outFilterMismatchNoverReadLmax 0.05 \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--twopassMode Basic \
--quantMode GeneCounts"

#########################

### EXECUTE COMMANDS ###
echo "### PHASE 1: STAGE-SPECIFIC FASTQ FILE ASSEMBLY AND ADAPTER TRIMMING ###"

# Gunzip any gzipped FASTQ files and move the files to the temp directory
fastq_file=$(find $sample_dir -name "$SAMPLE_NAME.*f*q*") 

if [[ $fastq_file = *.gz ]]
then
    echo "Unzipping fastq file $fastq_file to the sample temp directory"
    gunzip -c $fastq_file > $sample_dir/"$SAMPLE_NAME".fastq
fi

echo "Trimming transposase adapters with cutadapt"

cd $sample_dir

# Trim 0, 1, or 2 provided adapter sequences

keep_untrimmed="true"
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

# sizemarker_command="cutadapt \
    # -a 'AGCGTGTAGGGATCCAAA' \
    # -o "$SAMPLE_NAME"_trimmed.fastq \
    # "$SAMPLE_NAME"_adaptertrim.fastq"

# echo $sizemarker_command
# eval $sizemarker_command
cat "$SAMPLE_NAME"_adaptertrim.fastq > "$SAMPLE_NAME"_trimmed.fastq
python $python_dir/fastq_drop_short_reads.py . "$SAMPLE_NAME"_trimmed.fastq "$SAMPLE_NAME"_cleaned.fastq 16

echo "Removing trim intermediates"
rm -f "$SAMPLE_NAME"_trim1.fastq "$SAMPLE_NAME"_adaptertrim.fastq "$SAMPLE_NAME"_trimmed.fastq
echo "Adapter trimming complete."

echo "Filtering out low-complexity reads..."
python $python_dir/fastq_complexity_filter.py $sample_dir "$SAMPLE_NAME"_cleaned.fastq "$SAMPLE_NAME"_filtered.fastq 0.15
rm -f "$SAMPLE_NAME"_cleaned.fastq

echo "### PHASE 2: MAPPING READS TO THE GENOME WITH STAR ###"
mkdir -p $sample_dir/star
echo "STAR $star_params_allreads >& $SAMPLE_NAME.star.log"
eval "STAR $star_params_allreads >& $SAMPLE_NAME.star.log"

if [[ $fastq_file = *.gz ]]
then
   rm "$SAMPLE_NAME".fastq
fi

echo Generating bedgraph files...
samtools view -h star/Aligned.out.bam > "$SAMPLE_NAME".sam

python $python_dir/sam_calculate_coverages.py \
    -S "$SAMPLE_NAME".5p \
    -I "$SAMPLE_NAME".sam \
    -R TSS \
    -F $genome_fasta \
    --softclip_type 5p \
    --untemp_out G \
    --allow_naive \
    --minmatch 16


rm "$SAMPLE_NAME".sam 

bedtools sort -i "$SAMPLE_NAME".5p_TSS_plus.bedgraph > "$SAMPLE_NAME"_plus.5p.bedgraph
bedtools sort -i "$SAMPLE_NAME".5p_TSS_minus.bedgraph > "$SAMPLE_NAME"_minus.5p.bedgraph

bedtools sort -i "$SAMPLE_NAME".5p_TSS_plus_untemp.bedgraph > "$SAMPLE_NAME"_plus.uG.bedgraph
bedtools sort -i "$SAMPLE_NAME".5p_TSS_minus_untemp.bedgraph > "$SAMPLE_NAME"_minus.uG.bedgraph

rm "$SAMPLE_NAME".5p_TSS_plus.bedgraph \
    "$SAMPLE_NAME".5p_TSS_minus.bedgraph \
    "$SAMPLE_NAME".5p_TSS_plus_untemp.bedgraph \
    "$SAMPLE_NAME".5p_TSS_minus_untemp.bedgraph

rm "$SAMPLE_NAME"*coverage.bedgraph

# Calculate transcript_level coverage by parsing the genome coverage values
python $python_dir/bedgraph_genome_to_transcripts.py \
    --subset $annotation_subset \
    --output "$SAMPLE_NAME"_transcript.bedgraph \
    "$SAMPLE_NAME"_plus.5p.bedgraph \
    "$SAMPLE_NAME"_minus.5p.bedgraph \
    $annotation_gff \
    $genome_fasta

echo "Moving final files to results folder..."
mkdir -p $results_dir/EndMap/$SAMPLE_NAME
cp $sample_dir/*.bedgraph $results_dir/EndMap/$SAMPLE_NAME/
cat $sample_dir/comptable.tsv > $results_dir/EndMap/$SAMPLE_NAME/"$SAMPLE_NAME".comptable.tsv

echo Pipeline complete!

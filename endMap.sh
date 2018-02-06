#!/bin/bash

# EndMap
# Aligns RNA-seq reads from 5P or BODY experiments to a reference genome
# 1) Trims adapter sequences with cutadapt
# 2) Performs gapped alignment with STAR
# 3) Converts mapped reads to coverage values with readmapIO

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
if [ -z "$BIAS" ]
then
    BIAS=false
fi
BIAS=true

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
fastq_dir=${input_array[1]}    # Directory containing the fastq file
input_fastq=${input_array[2]}  # FASTQ filename(s), comma-separated
sample_name=${input_array[3]}  # Name of the sample
library_type=${input_array[4]} # Options: BODY, 5P. Type of the FASTQ file
read_type=${input_array[5]}    # Options: SE, PE, for single end or paired end libraries
adapter_str=${input_array[6]}  # comma-separated list of sequences to trim from the 3' end of reads
RAM=30
minimum_readlength=16
minimum_pairlength=16

echo "Parsed reference table:"
echo "Line number: $line_number"
echo "FASTQ directory: $fastq_dir"
echo "FASTQ file name: $input_fastq"
echo "Sample name: $sample_name"
echo "Library type: $library_type"
echo "Read type: $read_type"
echo "Adapter sequence(s): $adapter_str"

mkdir -p $temp_dir/star
temp_dir_s=$temp_dir/star/$sample_name
rm -rf $temp_dir_s
OPTIONS_star_global=$(cat $resource_dir/OPTIONS_star_global)

sample_dir=$temp_dir/$sample_name
rm -rf $sample_dir
mkdir -p $sample_dir

if [ $read_type == "SE" ]
then
    if [ $library_type == "5P" ]
    then
        star_params_allreads="$OPTIONS_star_global \
        --alignEndsType Local \
        --outFilterMatchNminOverLread 0.9 \
        --readFilesIn $sample_dir/"$sample_name"_cleaned.1.fastq \
        --outFileNamePrefix $sample_dir/star/"
    else
        star_params_allreads="$OPTIONS_star_global \
        --alignEndsType EndToEnd \
        --outFileNamePrefix $sample_dir/star/ \
        --readFilesIn $sample_dir/"$sample_name"_cleaned.1.fastq"
    fi
elif [ $read_type == "PE" ]
then
    star_params_allreads="$OPTIONS_star_global \
    --alignEndsType EndToEnd \
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

if [[ $library_type == "5P" ]]
then
    adapters=( $(echo $adapter_str | tr "," " ") )
    number_of_adapters=${#adapters[@]}
    echo "Adapters: ${adapters[@]}"
    if [ -z $adapters ]
    then
        echo "No adapter sequences specified. Skipping adapter trimming."
        cat "$sample_name".1.fastq > "$sample_name"_trimmed.fastq
    elif [[ $number_of_adapters -eq 2 ]]
    then
        first_trim_command="cutadapt \
            -a ${adapters[0]} \
            -o "$sample_name"_trim1.fastq \
            --untrimmed-output "$sample_name"_untrimmed1.fastq \
            "$sample_name".1.fastq"

        second_trim_command="cutadapt \
            -a ${adapters[1]} \
            -o "$sample_name"_trim2.fastq \
            "$sample_name"_untrimmed1.fastq"

        echo $first_trim_command
        eval $first_trim_command
        echo $second_trim_command
        eval $second_trim_command

        cat "$sample_name"_trim1.fastq "$sample_name"_trim2.fastq > "$sample_name"_trimmed.fastq
        rm -f "$sample_name"_trim1.fastq "$sample_name"_trim2.fastq "$sample_name"_untrimmed1.fastq

    elif [[ $number_of_adapters -eq 1 ]]
    then
        trim_command="cutadapt \
            -a ${adapters[0]} \
            -o "$sample_name"_adaptertrim.fastq \
            --untrimmed-output "$sample_name"_untrimmed1.fastq \
            "$sample_name".1.fastq"
            
        echo $trim_command
        eval $trim_command
        if [[ $keep_untrimmed == "true" ]]
        then
            cat "$sample_name"_adaptertrim.fastq "$sample_name"_untrimmed1.fastq > "$sample_name"_trimmed.fastq
        else
            cat "$sample_name"_adaptertrim.fastq > "$sample_name"_trimmed.fastq
        fi
    fi
      
    python $python_dir/fastq_drop_short_reads.py $sample_dir "$sample_name"_trimmed.fastq "$sample_name"_cleaned.fastq $minimum_readlength
    
    echo "Removing trim intermediates"
    rm -f "$sample_name"_trim1.fastq "$sample_name"_adaptertrim.fastq "$sample_name"_trimmed.fastq
    echo "Adapter trimming complete."

    echo "Filtering out low-complexity reads..."
    python $python_dir/fastq_complexity_filter.py $sample_dir "$sample_name"_cleaned.fastq "$sample_name"_cleaned.1.fastq 0.15
    rm -f "$sample_name"_cleaned.fastq
else
    #TODO: Improve adapter trimming behavior for BODY reads
    adapters=( $(echo $adapter_str | tr "," " ") )
    echo "Adapters: ${adapters[@]}"
    TN5_1=${adapters[0]}
    TN5_2=${adapters[1]}
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
        # directory,mate1,mate2,out1,out2,minlen=sys.argv[1:7]
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
echo "STAR $star_params_allreads >& $sample_name.star.log"
eval "STAR $star_params_allreads >& $sample_name.star.log"

# rm "$sample_name".fastq

echo Generating bedgraph files...
samtools view -H star/Aligned.out.bam | grep '@SQ' | sed 's/^@SQ\tSN://' | sed 's/\tLN:/\t/' > length.table

if [[ $BIAS == "true" ]]
then
    samtools view -h star/Aligned.out.bam > full.unsorted.sam
    python $python_dir/sam_subset.py \
        -F $genome_fasta \
        -A $annotation_gff \
        -I full.unsorted.sam \
        --feature CDS \
        --unique > sub.unsorted.sam
    samtools sort -O BAM -o sub.sorted.bam sub.unsorted.sam
    samtools index sub.sorted.bam
    samtools sort -O BAM -o full.sorted.bam full.unsorted.sam
    samtools index full.sorted.bam
    ##### Adapted from seqbias_pipe.sh ######
    # Wang et al. (2017) BMC Bioinformatics
    # https://github.com/txje/sequence-bias-adjustment
    
    # setup environment
    bias_pipe=$root_dir/scripts/sequence-bias-adjustment
    bias_outdir=$sample_dir/seqbias
    k=4
    mkdir -p $bias_outdir
    bam=sub.sorted.bam
    fullbam=full.sorted.bam
    bamnpy=$bias_outdir/$sample_name.bam.npy
    fullbamnpy=$bias_outdir/$sample_name.fullbam.npy
    baseline=$bias_outdir/$sample_name.read_50_baseline.csv
    kmer_bias=$bias_outdir/$sample_name.${k}mer_frequencies.csv
    tile_cov=$bias_outdir/$sample_name.tile_covariance.npy
    corrected_weights=$bias_outdir/$sample_name.${k}mer_adjusted.read_weights.npy
    corrected_full_weights=$bias_outdir/$sample_name.${k}mer_adjusted.read_weights.full.npy
    chroms=length.table
    read_len=$(samtools view $bam | head -n 1 | awk '{print $10}' - | wc | awk '{print $3}' -)
    if [ $read_len == 2 ]
    then
      read_len=20
      echo "$sample_name read length adjusted to: $read_len"
    fi
    # convert bam to npy array
    python $bias_pipe/bam2npy.py $bam $chroms $bamnpy
    python $bias_pipe/bam2npy.py $fullbam $chroms $fullbamnpy
    # compute nucleotide bias
    python $bias_pipe/compute_bias.py $bamnpy $genome_fasta $chroms 1 $bias_outdir/$sample_name.allele_frequencies.csv --read_len $read_len
    # compute k-mer baseline
    python $bias_pipe/compute_baseline.py $bamnpy $genome_fasta $chroms $k $baseline
    # compute k-mer bias
    python $bias_pipe/compute_bias.py $bamnpy $genome_fasta $chroms $k $kmer_bias --read_len $read_len
    # compute tile covariance matrix
    python $bias_pipe/correlate_bias.py $bamnpy $genome_fasta $chroms $kmer_bias $tile_cov
    # correct bias in subset
    python $bias_pipe/correct_bias.py $bamnpy $genome_fasta $chroms $baseline $kmer_bias $bias_outdir/$sample_name.${k}mer_adjusted.allele_frequencies.csv $corrected_weights $tile_cov --read_len $read_len
    # compute corrected nucleotide bias
    python $bias_pipe/compute_bias.py $corrected_weights $genome_fasta $chroms 1 $bias_outdir/$sample_name.${k}mer_adjusted.allele_frequencies.csv --read_len $read_len
    # comput corrected k-mer bias
    python $bias_pipe/compute_bias.py $corrected_weights $genome_fasta $chroms $k $bias_outdir/$sample_name.${k}mer_adjusted.${k}mer_frequencies.csv --read_len $read_len
    # apply bias correction to full bam file
    python $bias_pipe/correct_bias.py $fullbamnpy $genome_fasta $chroms $baseline $kmer_bias $bias_outdir/$sample_name.full.${k}mer_adjusted.allele_frequencies.csv $corrected_full_weights $tile_cov --read_len $read_len
    # output reweighted bam file (weights stored to XW tag)
    python $bias_pipe/npy2bam.py $chroms $corrected_full_weights $fullbam full.adjusted.bam --tag
    # resort adjusted bam file by read name
    samtools view -h -n full.adjusted.bam > $sample_name.sam
    # cleanup temp files
    rm -f sub.sorted.bam sub.unsorted.sam full.unsorted.sam $bamnpy $fullbamnpy
else
    samtools view -h star/Aligned.out.bam > $sample_name.sam
fi

if [[ $library_type == "5P" ]]
then
    python $python_dir/readmapIO.py \
        -S "$sample_name" \
        -I "$sample_name".sam \
        -R $library_type \
        -F $genome_fasta \
        --softclip_type 5p \
        --untemp_out A C G T \
        --secondary \
        --allow_naive
else
    python $python_dir/readmapIO.py \
        -S "$sample_name" \
        -I "$sample_name".sam \
        -R $library_type \
        -F $genome_fasta \
        --secondary \
        --allow_nonstranded \
        --allow_naive
fi

rm "$sample_name".sam 

bedtools sort -i "$sample_name"_"$library_type"_plus.bedgraph > "$sample_name"_plus.bedgraph
bedtools sort -i "$sample_name"_"$library_type"_minus.bedgraph > "$sample_name"_minus.bedgraph

if [[ $library_type == "5P" ]]
then
    bedtools sort -i "$sample_name"_"$library_type"_plus_untemp.bedgraph | \
        awk '{ printf $1"\t"$2"\t"$3"\t"$6"\n" }' > "$sample_name"_plus.uG.bedgraph
    bedtools sort -i "$sample_name"_"$library_type"_minus_untemp.bedgraph | \
        awk '{ printf $1"\t"$2"\t"$3"\t"$6"\n" }' > "$sample_name"_minus.uG.bedgraph
    # Calculate transcript_level coverage by parsing the genome coverage values
    python $python_dir/bedgraph_genome_to_transcripts.py \
        --subset $annotation_subset \
        --output "$sample_name"_transcript.bedgraph \
        "$sample_name"_plus.bedgraph \
        "$sample_name"_minus.bedgraph \
        $annotation_gff \
        $genome_fasta
else
    # Calculate transcript_level coverage by parsing the genome coverage values
    python $python_dir/bedgraph_genome_to_transcripts.py \
        --subset $annotation_subset \
        --output "$sample_name"_transcript.bedgraph \
        "$sample_name"_minus.bedgraph \
        "$sample_name"_plus.bedgraph \
        $annotation_gff \
        $genome_fasta
    
fi

rm "$sample_name"_"$library_type"_plus.bedgraph \
    "$sample_name"_"$library_type"_minus.bedgraph

# rm "$sample_name"*coverage.bedgraph

echo "Moving final files to results folder..."
mkdir -p $results_dir/EndMap/$sample_name
cp $sample_dir/*.bedgraph $results_dir/EndMap/$sample_name/
cat $sample_dir/comptable.tsv > $results_dir/EndMap/$sample_name/"$sample_name".comptable.tsv

echo Pipeline complete!

#!/bin/bash

function usage() {
cat <<HELP

### EndMap ###
Usage: ./endMap.sh [options] -L|--line <line_number>

Aligns RNA-seq reads from 5P and/or BODY experiments to a reference genome.
 1) Trims adapter sequences with cutadapt
 2) Performs gapped alignment with STAR
 3) Converts mapped reads to coverage values with readmapIO

Outputs a BEDGRAPH of mapped read 5' ends for both strands of the genome.

Optional arguments:
-R | --reference     Reference table (default: resources/reference.table)
-G | --genome        Genome FASTA file (default: resources/genome.fasta)
-A | --annotation    Transcript GFF file (default: resources/annotation.gff)
-B | --bias          Perform nucleotide bias correction (default: false)
--lmod               Load required modules with Lmod (default: false)
--ram                Amount of available RAM in gigabytes (default: 30)
--cpus               Number of cores available for multithreaded programs (default: 1)
--icomp              Minimum i-complexity score to filter reads before mapping (default: 0)

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

############################
# READING THE COMMAND LINE #
############################
# Taking the default environment above, add the commandline
# arguments (see read_cmdline.sh), and write a config file
if [ $# -eq 0 ]; then
    if [ -z "${PBS_ARRAY_INDEX}" ]
    then
        usage
        exit 1
    else
        JOB_NUMBER=${PBS_ARRAY_INDEX} # Imports the PBS job array number if it exists. Can be overridden with the commandline argument -J $JOB_NUMBER
    fi
else
    . $bash_dir/read_cmdline.sh
fi

# Set defaults for variables if they are not already in the environment from read_cmdline
if [ -z "$JOB_NUMBER" ]
then
    if [ -z "${PBS_ARRAY_INDEX}" ]
    then
        echo "ERROR: Please input a line number to process from reference.table"
        exit 1
    else
        JOB_NUMBER=${PBS_ARRAY_INDEX} # Imports the PBS job array number if it exists. Can be overridden with the commandline argument -J $JOB_NUMBER
    fi
fi
if [ -z "$GENOME_FASTA" ]
then
    GENOME_FASTA=$resource_dir/genome.fasta # If not already in environment, set as default value
fi
if [ -z "$REFERENCE_TABLE" ]
then
    REFERENCE_TABLE=$resource_dir/reference.table # If not already in environment, set as default value
fi
if [ -z "$ANNOTATION_GFF" ]
then
    ANNOTATION_GFF=$resource_dir/annotation.gff # If not already in environment, set as default value
fi
if [ -z "$LMOD" ]
then
    LMOD=0
fi
if [ -z "$CPUS" ]
then
    CPUS=1
fi
if [ -z "$RAM" ]
then
    RAM=30
fi
if [ -z "$BIAS" ]
then
    BIAS=0
fi
if [ -z "$ICOMP" ]
then
    ICOMP=0
fi

echo "Config settings:"
. $bash_dir/list_settings.sh

# Environment modules to load with Lmod (if option --lmod is passed)
REQUIRED_MODULES=( --bedtools --cutadapt --rna-star --samtools )
. $bash_dir/load_modules.sh
echo " "

input_mapper=$(sed -n "$JOB_NUMBER"p $REFERENCE_TABLE) #read mapping file
input_array=($input_mapper)

line_number=${input_array[0]}  # Line number of reference table (should match job number)
fastq_dir=${input_array[1]}    # Directory containing the fastq file
input_fastq=${input_array[2]}  # FASTQ filename(s), comma-separated
sample_name=${input_array[3]}  # Name of the sample
sample_type=${input_array[4]}  # Type of sample (samples of like type can be merged)
library_type=${input_array[5]} # Options: BODY, 5P. Type of the FASTQ file
read_type=${input_array[6]}    # Options: SE, PE, for single end or paired end libraries
adapter_str=${input_array[7]}  # comma-separated list of sequences to trim from the 3' end of reads
minimum_readlength=16
minimum_pairlength=16
sample_dir=$temp_dir/$sample_name.$library_type

echo "Parsed reference table:"
echo "Line number: $line_number"
echo "FASTQ directory: $fastq_dir"
echo "FASTQ file name: $input_fastq"
echo "Sample name: $sample_name"
echo "Sample type: $sample_type"
echo "Library type: $library_type"
echo "Read type: $read_type"
echo "Adapter sequence(s): $adapter_str"

mkdir -p $temp_dir/star
temp_dir_s=$temp_dir/star/$sample_name.$library_type
rm -rf $temp_dir_s
OPTIONS_star_global=$(cat $resource_dir/OPTIONS_star_global)

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
        echo "Multiple FASTQ files provided for single-end experiment. Merging..."
        fastq_split=$(echo $input_fastq | tr ',' ' ')
        fastq_array=($fastq_split)
        fastq_file_number=${#fastq_array[@]}
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
        gunzip -c $fastq_dir/$fastq_file > $sample_dir/"$sample_name"."$mate_number".fastq
    else
        echo "Copying fastq file $fastq_file to the sample temp directory"
        cat $fastq_dir/$fastq_file > $sample_dir/"$sample_name"."$mate_number".fastq
    fi
    if [[ $read_type == "SE"* ]] && [ $mate_number -gt 1 ]
    then
        cat $sample_dir/"$sample_name"."$mate_number".fastq >> $sample_dir/"$sample_name".1.fastq
        rm $sample_dir/"$sample_name"."$mate_number".fastq
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
    if [ $adapter_str == "nextera" ]
    then
        adapter_str="CTGTCTCTTATACACATCTGACGCTGCCGACGA,CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
    fi
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
      
    python $python_dir/fastq_drop_short_reads.py -O $sample_dir/"$sample_name"_cleaned.fastq --minlen $minimum_readlength $sample_dir/"$sample_name"_trimmed.fastq
    
    echo "Removing trim intermediates"
    rm -f "$sample_name"_trim1.fastq "$sample_name"_adaptertrim.fastq "$sample_name"_trimmed.fastq
    echo "Adapter trimming complete."
    
    if [ $(echo "$ICOMP > 0" | bc) -eq 1 ]
    then
        echo "Filtering out low-complexity reads..."
        python $python_dir/fastq_complexity_filter.py $sample_dir "$sample_name"_cleaned.fastq "$sample_name"_cleaned.1.fastq $ICOMP
    else
        cat "$sample_name"_cleaned.fastq > "$sample_name"_cleaned.1.fastq
    fi
    rm -f "$sample_name"_cleaned.fastq
else
    #TODO: Improve adapter trimming behavior for BODY reads
    if [ $adapter_str == "nextera" ]
    then
        adapter_str="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG,GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    fi
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
    ###########
    if [[ $read_type == "PE"* ]]
    then
        echo "cutadapt -a $TN5_1rc -o "$sample_name"_adaptertrim.2.fastq "$sample_name".2.fastq #rc(Tn5.1)"
        eval "cutadapt -a $TN5_1rc -o "$sample_name"_adaptertrim.2.fastq "$sample_name".2.fastq #rc(Tn5.1)"

        echo "Cleaning reads with fastq_drop_short_pairs.py"
        # directory,mate1,mate2,out1,out2,minlen=sys.argv[1:7]
        CMDline="python $python_dir/fastq_drop_short_pairs.py --minlen $minimum_pairlength -O $sample_dir/"$sample_name"_cleaned.1.fastq $sample_dir/"$sample_name"_cleaned.2.fastq $sample_dir/"$sample_name"_adaptertrim.1.fastq $sample_dir/"$sample_name"_adaptertrim.2.fastq"
        echo $CMDline
        eval $CMDline
        rm "$sample_name"_adaptertrim.1.fastq "$sample_name"_adaptertrim.2.fastq "$sample_name".1.fastq "$sample_name".2.fastq
    else
        echo "Cleaning reads with fastq_drop_short_reads.py"
        CMDline="python $python_dir/fastq_drop_short_reads.py -O $sample_dir/"$sample_name"_cleaned.1.fastq --minlen $minimum_readlength $sample_dir/"$sample_name"_adaptertrim.1.fastq"
        echo $CMDline
        eval $CMDline
        rm "$sample_name"_adaptertrim.1.fastq "$sample_name"_adaptertrim.2.fastq "$sample_name".1.fastq "$sample_name".2.fastq
    fi
    echo "Adapter trimming complete."
fi

echo "### PHASE 2: MAPPING READS TO THE GENOME WITH STAR ###"
mkdir -p $sample_dir/star
cd $sample_dir
echo "STAR $star_params_allreads >& $sample_name.star.log"
eval "STAR $star_params_allreads >& $sample_name.star.log"

# rm "$sample_name".fastq

echo Generating bedgraph files...
samtools view -H star/Aligned.out.bam | grep '@SQ' | sed 's/^@SQ\tSN://' | sed 's/\tLN:/\t/' > length.table

if [ $BIAS -eq 1 ]
then
    cd $sample_dir
    #TODO: Reconcile Python dependencies (2.7, 3.5)
    REQUIRED_MODULES=( --foss --python )
    . $bash_dir/load_modules.sh

    samtools view -h star/Aligned.out.bam > full.unsorted.sam
    python $python_dir/sam_subset.py \
        -F $GENOME_FASTA \
        -A $ANNOTATION_GFF \
        -I full.unsorted.sam \
        --feature CDS \
        --nucfreqs \
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
    k=3
    MARGIN=25
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
    
    REQUIRED_MODULES=( --pysam --numpy )
    . $bash_dir/load_modules.sh
    
    # convert bam to npy array
    python $bias_pipe/bam2npy.py $bam $chroms $bamnpy
    python $bias_pipe/bam2npy.py $fullbam $chroms $fullbamnpy
    # compute nucleotide bias
    python $bias_pipe/compute_bias.py $bamnpy $GENOME_FASTA $chroms 1 $bias_outdir/$sample_name.nucleotide_frequencies.csv --read_len $read_len --max_at_pos 5 --margin $MARGIN
    # compute k-mer baseline
    python $bias_pipe/compute_baseline.py $bamnpy $GENOME_FASTA $chroms $k $baseline --window_max $MARGIN
    # compute k-mer bias
    python $bias_pipe/compute_bias.py $bamnpy $GENOME_FASTA $chroms $k $kmer_bias --read_len $read_len --max_at_pos 5 --margin $MARGIN
    # compute tile covariance matrix
    python $bias_pipe/correlate_bias.py $bamnpy $GENOME_FASTA $chroms $kmer_bias $tile_cov --margin $MARGIN
    # correct bias in subset
    python $bias_pipe/correct_bias.py $bamnpy $GENOME_FASTA $chroms $baseline $kmer_bias $bias_outdir/unused_slot.csv $corrected_weights $tile_cov --read_len $read_len --max_at_pos 5 --margin $MARGIN
    # compute corrected nucleotide bias
    python $bias_pipe/compute_bias.py $corrected_weights $GENOME_FASTA $chroms 1 $bias_outdir/$sample_name.${k}mer_adjusted.nucleotide_frequencies.csv --read_len $read_len --max_at_pos 5 --margin $MARGIN
    # comput corrected k-mer bias
    python $bias_pipe/compute_bias.py $corrected_weights $GENOME_FASTA $chroms $k $bias_outdir/$sample_name.${k}mer_adjusted.${k}mer_frequencies.csv --read_len $read_len --max_at_pos 5 --margin $MARGIN
    # apply bias correction to full bam file
    python $bias_pipe/correct_bias.py $fullbamnpy $GENOME_FASTA $chroms $baseline $kmer_bias $bias_outdir/unused_slot.full.csv $corrected_full_weights $tile_cov --read_len $read_len --margin $MARGIN
    # output reweighted bam file (weights stored to XW tag)
    python $bias_pipe/npy2bam.py $chroms $corrected_full_weights $fullbam full.adjusted.bam --tag
    # resort adjusted bam file by read name
    
    REQUIRED_MODULES=( --samtools )
    . $bash_dir/load_modules.sh
    
    samtools sort -n full.adjusted.bam > $sample_name.sort.bam
    samtools view -h $sample_name.sort.bam > $sample_name.sam
    # cleanup temp files
    rm -f sub.unsorted.sam full.unsorted.sam $bamnpy $fullbamnpy
else
    samtools view -h star/Aligned.out.bam > $sample_name.sam
fi

REQUIRED_MODULES=( --foss --python )
. $bash_dir/load_modules.sh

if [ $BIAS -eq 1 ]
then
    if [[ $library_type == "5P" ]]
    then
        python $python_dir/readmapIO.py \
            -S "$sample_name" \
            -I "$sample_name".sam \
            -R $library_type \
            -F $GENOME_FASTA \
            --softclip_type 5p \
            --untemp_out G \
            --secondary \
            --allow_naive \
            --weighted
    else
        python $python_dir/readmapIO.py \
            -S "$sample_name" \
            -I "$sample_name".sam \
            -R $library_type \
            -F $GENOME_FASTA \
            --secondary \
            --allow_nonstranded \
            --allow_naive \
            --weighted    
    fi
else
    if [[ $library_type == "5P" ]]
    then
        python $python_dir/readmapIO.py \
            -S "$sample_name" \
            -I "$sample_name".sam \
            -R $library_type \
            -F $GENOME_FASTA \
            --softclip_type 5p \
            --untemp_out G \
            --secondary \
            --allow_naive
    else
        python $python_dir/readmapIO.py \
            -S "$sample_name" \
            -I "$sample_name".sam \
            -R $library_type \
            -F $GENOME_FASTA \
            --secondary \
            --allow_nonstranded \
            --allow_naive
    fi
fi

# rm "$sample_name".sam 

bedtools sort -i "$sample_name"_"$library_type"_plus.bedgraph > $sample_name.$library_type.plus.bedgraph
bedtools sort -i "$sample_name"_"$library_type"_minus.bedgraph > $sample_name.$library_type.minus.bedgraph

if [[ $library_type == "5P" ]]
then
    bedtools sort -i "$sample_name"_"$library_type"_plus_untemp.bedgraph | \
        awk '{ printf $1"\t"$2"\t"$3"\t"$4"\n" }' > $sample_name.uG.plus.bedgraph
    bedtools sort -i "$sample_name"_"$library_type"_minus_untemp.bedgraph | \
        awk '{ printf $1"\t"$2"\t"$3"\t"$4"\n" }' > $sample_name.uG.minus.bedgraph
else
    python $python_dir/bedgraph_combine.py \
        -i "$sample_name"_"$library_type"_plus_coverage.bedgraph "$sample_name"_"$library_type"_minus_coverage.bedgraph \
        -o $sample_name.$library_type.cov.bedgraph
fi

rm "$sample_name"_"$library_type"_plus.bedgraph \
    "$sample_name"_"$library_type"_minus.bedgraph

rm "$sample_name"*coverage.bedgraph

echo "Moving final files to results folder..."
mkdir -p $results_dir/EndMap/$sample_name
cp $sample_dir/*.bedgraph $results_dir/EndMap/$sample_name/

if [ $(echo "$ICOMP > 0" | bc) -eq 1 ]
then
    cat $sample_dir/comptable.tsv > $results_dir/EndMap/$sample_name/"$sample_name".comptable.tsv
fi

echo Pipeline complete!

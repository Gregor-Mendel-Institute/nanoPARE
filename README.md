# nanoPARE Analysis Tools
A data analysis pipeline companion to Schon, Kellner, et al.:  
["nanoPARE: Parallel analysis of RNA 5′ ends from low-input RNA"](https://genome.cshlp.org/content/28/12/1931)
  
<img width="200" alt="nanoPARE logo" src="/resources/images/nanoPARE_logo.png">
  
This repository contains a collection of tools for analyzing RNA 5′ end data from nanoPARE sequencing libraries. The method is designed as an extension of the widely used single-cell sequencing method [Smart-seq2](https://doi.org/10.1038/nmeth.2639) and allows the user to capture full-length transcript information afforded by Smart-seq2 in parallel with a library of 5' ends. By integrating the two pieces of information, the tools in this pipeline attempt to isolate the most relevant 5' ends in a genomic region by comparing the 5' end signal (5P reads) strength against a background model of full-length transcript data (BODY reads). The 5P features identified in the genome are then classified as capped transcription start sites or noncapped 5' ends. Finally, small RNA guided cleavage events are identified in the noncapped portion of the nanoPARE library. A summary of the included tools and their functions:  

  -**EndMap**: Align 5P and BODY FASTQ files to a reference genome.  
  -**EndGraph**: Identify 5P features by subtractive kernel density estimation.  
  -**EndClass**: Classify 5P features as capped or noncapped; label each feature according to a reference transcriptome.  
  -**EndMask**: Mask genomic regions with capped features and convert from genome coordinates to transcriptome coordinates.  
  -**EndCut**: Search for evidence of small RNA mediated cleavage in transcript-mapping noncapped 5P reads.  
  
These tools were written for Linux systems and optimized for a high-performance computing environment. Read alignment with EndMap uses the [STAR RNA-seq aligner](https://doi.org/10.1093/bioinformatics/bts635) and should be run with at least 30GB of RAM.  
  
Software requirements:  
  *STAR aligner 2.5+  
  *Python 3.6+  
  *Samtools 1.3+  
  *Bedtools 2.30+  
  *Cutadapt 1.9  
  
To download the repository, go to the desired destination folder and run:
```
git clone https://github.com/Gregor-Mendel-Institute/nanoPARE  
```
  
This will install the shell scripts described above, all default configuration files and Python utilities written for the pipeline, and a test dataset representing 1 megabase of the *Arabidopsis thaliana* reference genome and a subset of the data from our manuscript (GEO accession [GSE112869](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112869)). If your environment uses the module system [Lmod](https://github.com/TACC/Lmod), you should specify which modules to load for each piece of software above by modifying the file **resources/OPTIONS_lmod**. The steps below walk through the example data to show how each stage of the pipeline is run and give some information on the inputs and outputs of each step.  
  
  
### Setup (nanoPARE_setup.sh)  
**Usage: ./nanoPARE_setup.sh [options]**  
Options:  
  -G|--genome      Reference genome file in FASTA format
  -A|--annotation  Reference annotation file in GTF/GFF format
  -R|--reference   Reference table for input samples (see below)

Sets up a collection of files required to run the nanoPARE tools:  
  1) Writes a length.table recording the length of each chromosome in the genome FASTA  
  2) Generates a STAR genome index using the reference genome and transcriptome  
  3) Writes a BED file of putative TSO strand invasion sites based on mask_sequences.table  
  4) Writes a transcriptome FASTA file by combinding -G and -A  
  5) Parses reference annotations to identify 5'-terminal exons and single-exon transcripts  
  
This script must be run before performing any of the analysis steps for the first time. During setup, a few reference files are generated for the pipeline to recognize certain features in the reference genome (a multi-FASTA file), and reference transcriptome (a GTF or GFF3 formatted file that is indexed against the reference genome). *nanoPARE_setup.sh* uses these files to (1) find sites of potential strand invasion artifacts and (2) parse out a collection of 5'-most exons from the transcript annotations for comparison. If no files are provided, the setup will use **resources/genome.fasta** and **resources/annotation.gff**, which are part of the test dataset.  
  
To complete setup, you will also need to write an 8-column reference table that gives the pipeline all relevant information about the sequencing files you want to process. You can use the reference table in /resources/reference.table as a guide:  
  
| row number | directory | FASTQ filename | sample name | sample type | library type (5P or BODY) | sequencing run (SE or PE) | Adapter sequences |  
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |  
| 1 | resources/FASTQ | test5P_1.fq.gz | flower_1 | flower | 5P | SE | nextera |  
| 2 | resources/FASTQ | testBODY_1.fq.gz | flower_1 | flower | BODY | SE | nextera |  
| 3 | resources/FASTQ | test5P_2.fq.gz | flower_2 | flower | 5P | SE | nextera |  
| 4 | resources/FASTQ | testBODY_2.fq.gz | flower_2 | flower | BODY | SE | nextera |  
| 5 | resources/FASTQ | test5P_3.fq.gz | flower_3 | flower | 5P | SE | nextera |  
| 6 | resources/FASTQ | testBODY_3.fq.gz | flower_3 | flower | BODY | SE | nextera |  
  
All 5 tools will use this reference file as a lookup table. You can either modify the table in this repository (nanoPARE/resources/reference.table) or provide a filepath to your own table using the argument **-R|--reference** for any of the tools listed below. Each step was written so that samples can be processed in parallel, but downstream steps of the pipeline all use files generated in previous steps, so you will need to run each of them in order.  
  
### 1: EndMap (endMap.sh)  
**Usage: ./endMap.sh [options] -L|--line <line_number>**  
  
Aligns RNA-seq reads from 5P and/or BODY experiments to a reference genome.  
  1) Trims adapter sequences with cutadapt  
  2) Performs gapped alignment with STAR  
  3) Converts mapped reads to coverage values with readmapIO  
  
Outputs a BEDGRAPH of mapped read 5' ends for both strands of the genome.  
  
Optional arguments:
  --lmod   Load required modules with Lmod (default: false)
  --bias   Perform nucleotide bias correction (default: false)
  --ram    Amount of available RAM in gigabytes (default: 30)
  --cpus   Number of cores available for multithreaded programs (default: 1)
  --icomp  Minimum i-complexity score to filter reads before mapping (default: 0)
  
### 2: EndGraph (endGraph.sh)  
**Usage: ./endGraph.sh [options] -N|--name <sample_name>**  
  
Identification of end features through subtractive kernel density estimation  
  1) Determines a scaling factor to adjust read depths of 5P and BODY libraries  
  2) Smooths signal by fitting a Laplace kernel to END - BODY read values  
  3) Converts continuous regions of positive signal to features in a BED file  
  
Expects EndMap to be run and the results for sample_name in the directory /results/EndMap/sample_name  
  
Optional arguments:  
  --lmod       Load required modules with Lmod (default: false)  
  --cpus       Number of cores available for multithreaded programs (default: 1)  
  --rpm        Minimum RPM required to keep a feature (default: 5 for test data). Recommend 0.5 RPM for real libraries.  
  --kernel     Type of kernel to use for density estimation (default: laplace. options: gaussian, laplace)  
  --bandwidth  Bandwidth of kernel to use, in nucleotides (default: 15)  
  --fraglen    Mean fragment length of cDNA library, in nucleotides (default: 200)  
  
### 3: EndClass (endClass.sh)  
**Usage: ./endClass.sh [options] -T|--type <sample_type>**  
  
Labeling and classification of 5'-end features as capped or noncapped  
  1) Labels all 5P features based on their relationship to the nearest transcriptome annotation  
  2) Filters 5P features for only those replicable in >=2 experiments  
  3) Counts the proportion of reads in each 5P feature with upstream untemplated G (uuG)  
  4) Splits 5P features into "capped" (>=10% uuG) and "noncapped" (<10% uuG)  
  5) Calculates the number of full-length and truncated reads mapping to each gene  
  
Expects EndMap and EndGraph to be run before and their results to be in directories  
*/results/EndMap/sample_name* and */results/EndGraph/sample_name* respectively,  
for all samples of the chosen sample type that appear in reference.table  
  
Optional arguments:  
  --lmod      Load required modules with Lmod (default: false)  
  --cpus      Number of cores available for multithreaded programs (default: 1)  
  --uug       Minimum proportion uuG required to classify a feature as capped (default: 0.1)  
  --upstream  Maximum distance upstream (in nucleotides) to associate a 5P feature (default: 500)  
  
### 4: EndMask (endMask.sh)  
**Usage: ./endMask.sh [options] -T|--type <sample_type>**  
  
Masks reads belonging to capped features in the genome and prepares transcript bedgraph files for EndCut.  
  1) Identifies the transcript isoform with the most reads from each gene  
  2) Sets values for all 5P reads within a capped feature to 0  
  3) Writes a transcript-indexed bedgraph file of cap-masked reads with all dominant transcripts  
  
Expects EndMap, EndGraph, and EndClass to be run before and their results to be in directories  
*/results/EndMap/sample_name*, */results/EndGraph/sample_name*, and */results/EndClass/sample_type*  
respectively, for all samples of the chosen sample_type that appear in reference.table  
  
Optional arguments:  
  --lmod  Load required modules with Lmod (default: false)  
  --cpus  Number of cores available for multithreaded programs (default: 1)  
  --mask  Alternative sample_type to use for cap masking (default: sample_type)  
  
### 5: EndCut  
Sequences from miRNAs and tasiRNAs annotated in TAIR10 or miRBase21 (Lamesch et al. 2012; Kozomara and Griffiths-Jones 2013) were selected (i.e. anno.mir.tas.fa) and randomized one thousand times each by the python script sRNA_shuffler.py to produce anno.mir.tas.i.fa files; where i is an integer between 0 and 999). For annotated miRNA, tasiRNA and the corresponding 1,000 randomized variants for each miRNA/tasiRNA, GSTAr.pl (https://github.com/MikeAxtell/GSTAr) was used to predict target sites in transcript models annotated as protein-coding genes, transposable element genes or other RNAs (i.e. TAIR10_pc_teg_other_transcripts.fasta). Target sites were determined based on the level of complementarity between sRNAs and transcripts computed using previously developed criteria based on the frequency and position of the miRNA-target duplex mismatches (i.e. Allen scores) (Allen et al. 2005) As described above, nanoPARE data was processed by EndMask to exclude capped regions of transcripts from further analyses. Publicly available PARE datasets were downloaded from the Sequence Read Archive (NCBI) (Supplementary Data S1), but alignments overlapping capped features were not excluded from downstream analyses. 
Predicted target sites and EndGraph output were used by EndCut_step1.sh to quantify the number of reads at predicted target sites and in adjacent 20 nt or 50 nt regions on the sense strand of the same transcript. Adjacent sites within one nucleotide of predicted cleavage sites were not considered in order to not penalize sites for sRNA isoforms with slightly offset target recognition sites. The local enrichment of nanoPARE read 5′ ends at predicted cleavage sites relative to surrounding transcribed regions, or fold-changes, were calculated by dividing the numbers of nanoPARE read 5′ ends at predicted cleavage sites + 1 by the maximum numbers of reads in adjacent transcript regions + 1. Allen scores were also assigned to each predicted cleavage site detected. For each randomized sRNA control set, EndCut_step2.R computed empirical cumulative distribution functions of fold-changes (ECDFFC) and Allen scores (ECDFAS). These were then used as null models to test whether the observed cleavage site fold-changes were not equal to or lesser than ECDFFC, as well as if the observed site Allen scores were not equal to or greater than ECDFAS. Final P-values were computed for each site by combining these two P-values using Fisher’s combined probability test, and then adjusted for multiple testing using the Benjamini- and Hochberg method. For our analyses, we selected predicted cleavage sites with adjusted P-values < 0.05, fold-changes > 1.0 and that were also represented by at least one read per ten million transcriptome-mapping reads. 

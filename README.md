# nanoPARE Analysis Tools
A data analysis pipeline companion to Schon, Kellner, et al.:  
"nanoPARE: Parallel analysis of RNA 5′ ends from low-input RNA"  
  
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
  *Python 3.5+  
  *Samtools 1.3+  
  *Bedtools 2.26  
  *Cutadapt 1.9  
  
To download the repository, go to the desired destination folder and run:
```
git clone https://github.com/Gregor-Mendel-Institute/nanoPARE  
```
  
This will install the shell scripts described above, all default configuration files and Python utilities written for the pipeline, and a test dataset representing 1 megabase of the *Arabidopsis thaliana* reference genome and a subset of the data from our manuscript (GEO accession [GSE112869](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112869)). If your environment uses the module system [Lmod](https://github.com/TACC/Lmod), you should specify which modules to load for each piece of software above by modifying the file **resources/OPTIONS_lmod**. The steps below walk through the example data to show how each stage of the pipeline is run and give some information on the inputs and outputs of each step.  
  
  
### Setup (nanoPARE_setup.sh)  
INPUTS:  
  **-G|--genome** [file] (genome sequence in FASTA format)  
  **-A|--annotation** [file] (gene annotation file in GTF/GFF3 format)  
  
This script must be run before performing any of the analysis steps for the first time. During setup, a few reference files are generated for the pipeline to recognize certain features in the reference genome (a multi-FASTA file), and reference transcriptome (a GTF or GFF3 formatted file that is indexed against the reference genome). *nanoPARE_setup.sh* uses these files to (1) find sites of potential strand invasion artifacts and (2) parse out a collection of 5'-most exons from the transcript annotations for comparison.  
  
To complete setup, you will also need to write an 8-column reference table that gives the pipeline all relevant information about the sequencing files you want to process. You can use the reference table in /resources/reference.table as a guide:  
  
| row number | directory | FASTQ filename | sample name | sample type | library type (5P or BODY) | sequencing run (SE or PE) | Adapter sequences |  
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |  
| 1 | resources/FASTQ | test5P_1.fq | flower_1 | flower | 5P | SE | nextera |  
| 2 | resources/FASTQ | testBODY_1.fq | flower_1 | flower | BODY | SE | nextera |  
| 3 | resources/FASTQ | test5P_2.fq | flower_2 | flower | 5P | SE | nextera |  
| 4 | resources/FASTQ | testBODY_2.fq | flower_2 | flower | BODY | SE | nextera |  

  
All 5 tools will use this reference table as a lookup. You can either modify the table in this repository (nanoPARE/resources/reference.table) or provide a filepath to your own table using the argument **-R|--reference** for any of the tools listed below. Each step was written so that samples can be processed in parallel, but downstream steps of the pipeline all use files generated in previous steps, so you will need to run each of them in order.  
  
### 1: EndMap (endMap.sh)  
INPUTS:  
  **-L|--line** [number] Line number in the reference table of the sample you want to process  
OPTIONAL:  
  **-R|--reference** [file] A reference table as formatted above  
  **-I|--icomp** [float] minimum per-nucleotide i-complexity (default=0.15)
  
EndMap first trims the appropriate adapter sequences using *cutadapt*. By default these sequences are based on the Tn5 adapter sequences of the Nextera Tagmentation kit. To prevent reads with low sequence complexity from mapping to the genome, the I-complexity (Becher and Heiber 2012) of each FASTQ read is calculated by the Python script *fastq_complexity_filter.py*, and reads with a per-nucleotide I-complexity <0.15 are removed prior to mapping. The remaining reads are then aligned to the genome with *STAR*. The alignment settings are optimized for Arabopsis by default, you can change these settings in the file **resources/OPTIONS_star_global**. In addition to these global settings, mapping behavior differs slightly when BODY and 5P libraries. BODY reads are mapped using --alignEndsType EndToEnd, and softclipping is allowed at the ends of 5P reads using --alignEndsType Local and --outFilterMatchNminOverLread 0.9.  
After alignment, 5′ end read depth for each position in the genome is calculated from the BAM file using the Python script *readmapIO.py*. Each BAM file contains reads that mapped from 1-100 times to the reference genome. For reads that mapped to more than one location, the most likely locus of origin is inferred via a “rich-get-richer” algorithm similar to that employed by the software [MuMRescue](https://academic.oup.com/bioinformatics/article/25/19/2613/180391). ReadmapIO begins by calculating the coverage depth of uniquely mapping reads at each position in the genome. Multimappers are then binned by their mapping multiplicity (i.e. a read that maps to 10 locations in the genome has a multiplicity of 10). Beginning with a multiplicity of 2, all reads in that bin are sorted from lowest possible genomic position to highest, and each read is assigned in a multistep process: if at least one mapping position has at least one existing read, the read is considered “unambiguous”, and is assigned proportionally to its mapping locations using the formula <img alt="P_{i}=\frac{C_{i}}{\sum_{j=1}^{n} C_{j}}" src="/resources/images/readmapIO_eq1.gif">,  where *Pi* is the proportion of reads assigned to mapping location *i*, *Ci* is the total existing read coverage assigned to the genomic positions that comprise location *i*, and *n* is the number of mapping locations for the read. If the existing read coverage at all locations is 0, that read is not yet assigned. After examining all reads in the bin, the bin is sorted from highest genomic position to lowest and the process is repeated, until no more unambiguous reads can be identified. At this point, all remaining reads in the bin are assigned with equal weighting, or *Pi=1/n*. Bins are assigned this way in order from a multiplicity of 2 to 100. After all reads are assigned, readmapIO outputs a bedgraph file of assigned read 5′ end counts for the plus and minus strand of the genome. If the library type is 5P, readmapIO also outputs bedgraph files of all nucleotides softclipped from the 5′ end of reads. These will be referred to hereafter as upstream untemplated nucleotides (uuNs).

  
### 2: EndGraph (endGraph.sh)  
INPUTS:  
  **-N|--name** [string] sample name with both a 5P and BODY library to process from the reference table  
OPTIONAL:  
  **-R|--reference** A reference table as formatted above  
  
Discrete 5P features are identified genome-wide via subtractive kernel density estimation. Bedgraph files output from EndMap corresponding to a sample’s 5P and BODY libraries were evaluated together. First, strand invasion artifacts are masked based on the masking file that was generated during *nanoPARE_setup.sh*. Then, a scaling factor (S) is estimated to normalize the read depth of the 5P library against the BODY library using the formula:
<img alt="S=\frac{2F*10^{6}}{\sum_{i=1}^{n}TPM_i*L_i))}*\frac{R_B}{R_E}" src="/resources/images/endgraph_eq1.gif">,
where *n* is the total number of transcripts, *TPMi* is the abundance of a transcript in transcripts per million, *Li* is the length of a transcript in nucleotides, *F* is the mean fragment length of the BODY library, *RB* is the total number of mapped BODY reads, and *RE* is the total number of mapped 5P reads. Then, a Laplace kernel with a bandwidth of 15 nucleotides is fit over the set of values *(ER * S) – BR*, where *ER* is the set of 5P end read counts and *BR* is the set of body read counts. Regions of continuous positive density were extracted and written as discrete features to a bed file. The formula above attempts to normalize the two libraries by estimating the expected ratio of cDNA fragments  that contain 5’ ends based on the assumption that all signal is derived from full-length transcripts (as defined by the reference transcript annotations). By making this assumption, spurious 5P signal downstream of the true 5’ terminus should be weaker than it needs to be to outweigh the BODY signal in that location.  
  
### 3: EndClass (endClass.sh)  
INPUTS:  
  **-T|--type** [string] sample type with at least two samples included in the reference table  
OPTIONAL:  
  **-R|--reference** [file] A reference table as formatted above  
  **--min_uug** [float] Minimum percent upstream untemplated G to classify a feature as "capped" (default=0.1)
  
If a 5P experiment was designed with multiple replicates, EndClass merged all 5P features that could be reproducibly identified in ≥2 replicates. Then, the presence of a m7G cap was predicted for each replicable feature by calculating the proportion of reads containing upstream untemplated guanosine (uuG). A feature was considered capped if ≥10% of all reads from a sample type that map within the feature contained uuG, otherwise the feature was considered noncapped.  
  
### 4: EndMask (endMask.sh)  
INPUTS:  
  **-T|--type** [string] sample type with at least two samples included in the reference table  
OPTIONAL:  
  **-R|--reference** A reference table as formatted above  
  
EndMask prepared a bedgraph file of 5P read positions relative to the start site of the dominant isoform of each gene in the reference annotation. Dominant isoforms were defined as the transcript isoform containing the most mapped reads. For nanoPARE libraries, this transcript-level bedgraph was generated with a cap-masked input in which 5P reads contained within replicable capped 5P features were discarded.
  
### 5: EndCut (endCut.sh)  
INPUTS:  
  **-N|--name** [string] sample name to process from the reference table  
OPTIONAL:  
  **-R|--reference** A reference table as formatted above  
  
Sequences from miRNAs and tasiRNAs annotated in TAIR10 or miRBase21 (Lamesch et al. 2012; Kozomara and Griffiths-Jones 2013) were selected (i.e. anno.mir.tas.fa) and randomized one thousand times each by the python script sRNA_shuffler.py to produce anno.mir.tas.i.fa files; where i is an integer between 0 and 999). For annotated miRNA, tasiRNA and the corresponding 1,000 randomized variants for each miRNA/tasiRNA, GSTAr.pl (https://github.com/MikeAxtell/GSTAr) was used to predict target sites in transcript models annotated as protein-coding genes, transposable element genes or other RNAs (i.e. TAIR10_pc_teg_other_transcripts.fasta). Target sites were determined based on the level of complementarity between sRNAs and transcripts computed using previously developed criteria based on the frequency and position of the miRNA-target duplex mismatches (i.e. Allen scores) (Allen et al. 2005) As described above, nanoPARE data was processed by EndMask to exclude capped regions of transcripts from further analyses. Publicly available PARE datasets were downloaded from the Sequence Read Archive (NCBI) (Supplementary Data S1), but alignments overlapping capped features were not excluded from downstream analyses. 
Predicted target sites and EndGraph output were used by EndCut_step1.sh to quantify the number of reads at predicted target sites and in adjacent 20 nt or 50 nt regions on the sense strand of the same transcript. Adjacent sites within one nucleotide of predicted cleavage sites were not considered in order to not penalize sites for sRNA isoforms with slightly offset target recognition sites. The local enrichment of nanoPARE read 5′ ends at predicted cleavage sites relative to surrounding transcribed regions, or fold-changes, were calculated by dividing the numbers of nanoPARE read 5′ ends at predicted cleavage sites + 1 by the maximum numbers of reads in adjacent transcript regions + 1. Allen scores were also assigned to each predicted cleavage site detected. For each randomized sRNA control set, EndCut_step2.R computed empirical cumulative distribution functions of fold-changes (ECDFFC) and Allen scores (ECDFAS). These were then used as null models to test whether the observed cleavage site fold-changes were not equal to or lesser than ECDFFC, as well as if the observed site Allen scores were not equal to or greater than ECDFAS. Final P-values were computed for each site by combining these two P-values using Fisher’s combined probability test, and then adjusted for multiple testing using the Benjamini- and Hochberg method. For our analyses, we selected predicted cleavage sites with adjusted P-values < 0.05, fold-changes > 1.0 and that were also represented by at least one read per ten million transcriptome-mapping reads. 

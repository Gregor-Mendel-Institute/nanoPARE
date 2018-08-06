# nanoPARE
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
  
Additional software requirements:  
  -STAR aligner 2.5+  
  -Python 3.5+  
  -Samtools 1.3+  
  -Bedtools v2.17  
  
To download the repository, go to the desired destination folder and run:
```
git clone https://github.com/Gregor-Mendel-Institute/nanoPARE  
```
  
This will install the shell scripts described above, all default configuration files and Python utilities written for the pipeline, and a test dataset representing 1 megabase of the *Arabidopsis thaliana* reference genome and a subset of the data from our manuscript (GEO accession [GSE112869](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112869))

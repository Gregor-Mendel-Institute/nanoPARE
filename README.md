# EndGraph
Subtractive kernel density modeling of RNA 5' and/or 3' ends.

  A variety of methods have been developed for genome-wide detection of RNA 5' ends or 3' ends.
  Knowing the exact position of transcription initiation or polyadenylation is necessary to investigate functional DNA regulatory elements immediately upstream or downstream of a gene, as well as RNA elements in 5' or 3' UTRs. In addition to sites of transcription initiation, which in eukaryotes produce an RNA molecule containing a 5'-methylguanosine cap, uncapped 5' ends of RNA can also be informative, such as microRNA-targeted sites of endonucleolytic cleavage. Methods exist to capture 5' capped RNA, uncapped RNA, or capped and uncapped species together. However, reads in these libraries are often of limited value without the context of the full-length RNA.  
  This pipeline takes a pair of deep sequencing datasets from full-length RNA and from 5'/3' ends, then produces a subtractive model based on kernel density estimation. This model can be applied genome-wide to identify RNA end variations that are significant relative to the background expression of the gene, allowing for high-confidence annotation of novel extended or truncated ends, and for significant uncapped 5' ends in mixed datasets. EndGraph is optimized for use with a paired Smart-seq2 library and nanoPARE library from the same cDNA.
  
Software requirements:  
  -Bedtools  
  -Python 2.7+  
  -Samtools  

  
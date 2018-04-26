#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N get.dists
#PBS -J 1-12
#PBS -j oe
#PBS -o logs/get.dists.^array_index^.log
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=8gb

#set variables
#export NAME_SAMPLE=fb1_1to1
#export SAMPLE=fb.W
export outDir_s=transcript_bedgraph_capmasked
export outDir_a=W.transcript.capmasked
#########################################
export WSEQ_s=mRNA/180228/
export WSEQ=/lustre/scratch/users/michael.nodine/seq/$WSEQ_s
export outDir=$WSEQ/results/$outDir_s/
export dataRoot=WSEQ
#export annoRoot=/lustre/scratch/users/michael.nodine/annos/
export TEST=anno.mir.tas.fa.GSTAr
export SHUFF=anno.mir.tas
#export BG=$outDir_s/${NAME_SAMPLE}.$outDir_a.bedgraph
#export TF=$WSEQ/bedFiles/dominant_transcript_tables/${SAMPLE}.dominant_transcript_lengths.tsv
########################################

# build array index
case ${PBS_ARRAY_INDEX} in
        1) NAME='fb1_1to1';;
        2) NAME='fb2_1to1';;
        3) NAME='fb3_1to1';;
        4) NAME='xrn1_fb_1';;
        5) NAME='xrn1_fb_2';;
        6) NAME='xrn1_fb_3';;
        7) NAME='xrn4_M_1';;
        8) NAME='xrn4_M_2';;
        9) NAME='xrn4_M_3';;
        10) NAME='d234_fb4';;
        11) NAME='d234_fb5';;
        12) NAME='d234_fb6';;
esac

#mkdir $WSEQ/results/distances/$outDir_s/$NAME_SAMPLE/

#module load BEDTools/2.26.0-foss-2017a
module load BEDTools/v2.17.0-goolf-1.4.10

#Step 1 (single job): select GSTAr entries with Allen Scores between 0 and 6 and generate bedfile of slice site with Allen Score
awk 'NR > 8' $WSEQ/GSTAr/$TEST | awk -v OFS='\t' '{if ($9 >= 0 && $9 <= 6.0 && ($5 - 101 >0)) print $2,$5 - 101,$5 + 100,$1,$9; else if ($9 >= 0 && $9 <= 6.0) print $2,0,$5 + 100,$1,$9}' | sort -k1,1 -k2,2n > $WSEQ/results/distances/$outDir_s/$TEST.pred.sites.100nt.bed

#merge identical sites in bed; same site can regulated by multiple sRNAs/isoforms and should not be considered separately
bedtools merge -i $WSEQ/results/distances/$outDir_s/$TEST.pred.sites.100nt.bed -c 4,5 -o collapse,min > $WSEQ/results/distances/$outDir_s/$TEST.pred.sites.100nt.merged.bed

#Step 2 (parallel jobs):
#get overlapping sites in nanoPARE bedgraph file
awk 'BEGIN {print "Transcript\tsRNA\tAllen.Score\tPred.site\tDet.site\tDiff\tReads"}' > $WSEQ/results/distances/$outDir_s/${NAME}.pred.sites.100nt.overlap.tsv

bedtools intersect -wb -wa -a $WSEQ/results/distances/$TEST.pred.sites.100nt.bed -b $WSEQ/bedFiles/$outDir_s/${NAME}.$outDir_a.bedgraph > $WSEQ/results/distances/$outDir_s/${NAME}.pred.sites.100nt.bed.overlap 

awk -v OFS='\t' '{print $1,$4,$5,$3-101,$7,$7-($3-101),$9}' $WSEQ/results/distances/$outDir_s/${NAME}.pred.sites.100nt.bed.overlap >> $WSEQ/results/distances/$outDir_s/${NAME}.pred.sites.100nt.overlap.tsv


#remove intermediate files
rm $WSEQ/results/distances/$outDir_s/${NAME}.pred.sites.100nt.bed.overlap

echo 'nanoPARE_get.dists.sh complete!'

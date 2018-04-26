#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N bed.norm
#PBS -j oe
#PBS -o logs/bed.norm.log
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=2:mem=8gb

############set variables##############
export outDir_s=transcript_bedgraph_capmasked
export outDir_a=W.transcript.capmasked
export WSEQ_s=mRNA/180228/
export WSEQ=/lustre/scratch/users/michael.nodine/seq/$WSEQ_s
export NAME='fb1_1to1'
########################################

export BG=$WSEQ/bedFiles/$outDir_s/${NAME}.$outDir_a.bedgraph
export BG_norm=$BG.norm

echo 'Number of lines in original bedgraph:'
wc $BG

##get total number of transcriptome-mapping reads from BG file and use it to normalize BG hit-normalized reads
TMR=$(awk -v OFS='\t' '{sum+=$4}END{print sum}' $BG)
MTMR=$(awk '{print +$1 }' <<< $TMR)
MTMR=$(bc -l <<< $MTMR/1000000)

echo 'Total number of hit-normalized transcriptome-mapping reads (in millions:',$MTMR

awk -v OFS='\t' -v MTMR="$MTMR" '{print $1,$2,$3,$4,$4/MTMR}' $BG > $BG_norm

echo 'Number of lines in normalized bedgraph:'
wc $BG_norm

echo 'nanoPARE_EndCut.step1_v02_normalize.sh complete!'

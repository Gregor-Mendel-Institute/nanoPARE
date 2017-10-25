import sys
import os
import re
import math
import argparse
import fasta_utils as fu
import multiprocessing as mp
from collections import Counter

desc = (
    "Takes a transcript-level BEDGRAPH and outputs only the"
    "most covered isoform of each transcript."
    "Optionally outputs a table of transcript lengths for"
    "the transcripts included in the output bedgraph."
)
parser = argparse.ArgumentParser(description=desc)

# add arguments to the ArgumentParser
parser.add_argument(
    '-I', '--input', dest='input', type=str, 
    help='input bedgraph coverage file',
    required=True
)
parser.add_argument(
    '-F', '--fasta', dest='fasta', type=str, 
    help='input FASTA file',
    required=True
)
parser.add_argument(
    '-O', '--output', dest='output', type=str, 
    help='output table filepath',
    default='dominant_transcript.bedgraph'
)
parser.add_argument(
    '-L', '--lengths', type=str, default=None,
    help='filepath to output lengths table'
)
args = parser.parse_args()

def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]


def write_bedgraph_from_dict(input, output_filename, digits=2, parallel=False):
    """Writes unsorted BEDGRAPH to output_filename from input dict"""
    
    def writer(merge_filename,shared_queue,stop_token):
        """Initializes a writer to output the mp queue."""
        dest_file = open(merge_filename,'w')
        while True:
            line = shared_queue.get()
            if line == stop_token:
                dest_file.close()
                return
            dest_file.write(line)
    
    
    def generate_bedgraph_lines(values_dict, chrom, queue, parallel=parallel):
        """Converts dict of values to bedgraph format"""
        start = 0
        prevpos = 0
        prevcount = None
        chromlen = len(input[chrom])
        position = None
        # Iterate through every position with values in the dict
        for position in sorted(list(values_dict.keys())):
            count = round(values_dict[position],digits)
            if count != prevcount or int(position) > 1 + prevpos:
                # The newly encountered value is not a continuation
                # of the previous value. Write the old run and start another.
                if prevcount > 0:
                    line_to_write = '\t'.join(
                        [
                            str(i) for i in [
                                chrom,
                                start,
                                prevpos + 1,
                                prevcount
                            ]
                        ]
                    ) + '\n'
                    
                    if parallel:
                        # Write the old run to outfile
                        queue.put(line_to_write)
                    else:
                        # Write the old run to outfile
                        queue.write(line_to_write)
                        
                start = position
            prevcount = count
            prevpos = int(position)
        
        if position and prevcount > 0:
            line_to_write = '\t'.join(
                [
                    str(i) for i in [
                        chrom,
                        start,
                        prevpos + 1,
                        prevcount
                    ]
                ]
            ) + '\n'
            
            if parallel:
                queue.put(line_to_write)
            else:
                queue.write(line_to_write)
    
    
    if parallel:
        queue = mp.Queue()
        STOP_TOKEN = "FINISHED"
        writer_process = mp.Process(
            target=writer,
            args=(output_filename,queue,STOP_TOKEN)
        )
        writer_process.start()
        all_threads = []

        for chrom in sorted(list(input.keys())):
            all_threads.append(
                mp.Process(
                    target=generate_bedgraph_lines,
                    args=(
                        input[chrom],
                        chrom,
                        queue
                    )
                )
            )
        for i in range(len(all_threads)):
            all_threads[i].start()
        while len(mp.active_children()) > 1:
            time.sleep(1)
        queue.put("FINISHED")
        while len(mp.active_children()) > 0:
            time.sleep(1)
    else:
        queue = open(output_filename, 'w')
        for chrom in sorted(list(input.keys())):
            generate_bedgraph_lines(
                input[chrom],
                chrom,
                queue,
                parallel = False
            )
        queue.close()

#'chromosomes' contains the lengths of all chromosomes the that BEDGRAPH contains values for.

genome = fu.import_genome(args.fasta)
chromosomes = dict(
    [
        (ID,len(transcript))
        for ID,transcript in genome.items()
    ]
)

#'coverage' is a nested dictionary of float vectors for each nucleotide in the genome.
# Contains a dictionary for each BEDGRAPH file with values at each position.
coverage = {}
graph = open(args.input)
coverage = {}

for line in graph:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    if chrom not in coverage:
        coverage[chrom] = {}
    
    for i in range(int(start),int(end)):
        coverage[chrom][i] = count

# Generate a collection of genes and their isoforms
# by splitting each isoform name on the '.'
genes = {}
for i in genome.keys():
    gene_id = i.split('.')[0]
    genes[gene_id] = genes.get(gene_id,[]) + [i]

dominant_transcripts = []
for g in genes.keys():
    isos = sorted(genes[g])
    dom = None
    highest = 0
    for i in isos:
        v = 0
        if i in coverage:
            c = coverage[i].values()
            if len(c) == 0:
                v = 0
            else:
                v = sum(c)
        if v > highest:
            dom = i
    
    if dom:
        dominant_transcripts.append(dom)

# Generate a new dictionary of only dominant transcript coverages
dominant_coverage = dict(
    [
        (k,v)
        for k,v in coverage.items()
        if k in dominant_transcripts
    ]
)

# OUTPUT RESULTS
if args.lengths:
    length_table = open(args.lengths, 'w')
    for d in sorted(dominant_transcripts):
        length_table.write(
            '{}\t{}\n'.format(
                d,
                len(genome[d])
            )
        )
    length_table.close()
    
write_bedgraph_from_dict(
    dominant_coverage,
    args.output,
    digits=2,
    parallel=False
)


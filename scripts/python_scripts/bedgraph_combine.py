import sys
import os
import re
import math
import argparse
import multiprocessing
from collections import Counter

desc = (
    "Combines the values from two or more bedgraph files."
)
parser = argparse.ArgumentParser(description=desc)

# add arguments to the ArgumentParser
parser.add_argument(
    '-i', '--input', dest='INPUT', type=str, 
    help='input bedgraph filepath(s)',
    default=[], nargs='+'
)
parser.add_argument(
    '-s', '--scale', dest='SCALE', type=float, 
    help='(optional) A scaling weight for each input',
    default=[], nargs='+'
)
parser.add_argument(
    '-o', '--output', dest='OUTPUT', type=str, 
    help='output bedgraph filepath',
    default='out.bedgraph'
)

args = parser.parse_args()
if args.SCALE:
    if len(args.SCALE) != len(args.INPUT):
        print("ERROR: Number of scaling weights must match number of inputs.")
        sys.exit(1)

def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]

def write_bedgraph_from_dict(input, output_filename, digits=2, parallel=False, multi_key=False):
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
    
    
    def generate_bedgraph_lines(values_dict, chrom, queue, digits=digits, parallel=parallel, multi_key=multi_key):
        """Converts dict of values to bedgraph format"""
        start = 0
        prevpos = 0
        prevcount = None
        chromlen = len(input[chrom])
        position = None
        # Iterate through every position with values in the dict
        for position in sorted(list(values_dict.keys())):
            if multi_key:
                all_values = [str(round(values_dict[position].get(k,0),digits)) for k in multi_key]
                count = '\t'.join(all_values)
            else:
                count = round(values_dict[position],digits)
            
            if count != prevcount or int(position) > 1 + prevpos:
                # The newly encountered value is not a continuation
                # of the previous value. Write the old run and start another.
                if prevcount and prevcount != 0:
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
        
        if position and prevcount and prevcount != 0:
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
        
        if multi_key:
            first_line = '#chrom\tstart\tend\t{}\n'.format(
                '\t'.join([str(i) for i in multi_key])
            )
            queue.put(first_line)
        
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
        if multi_key:
            first_line = '#chrom\tstart\tend\t{}\n'.format(
                '\t'.join([str(i) for i in multi_key])
            )
            queue.write(first_line)
        
        for chrom in sorted(list(input.keys())):
            generate_bedgraph_lines(
                input[chrom],
                chrom,
                queue,
                parallel = False
            )
        queue.close()

    
coverage = {}

if args.INPUT:
    bedgraph_outfile=open(args.OUTPUT,'w')
    for filenumber in range(len(args.INPUT)):
        bedgraph_file=open(args.INPUT[filenumber])
        scale = float(1)
        if args.SCALE:
            scale = args.SCALE[filenumber]
        
        for line in bedgraph_file:
            chrom,start,end,count = line.rstrip().split()
            count = float(count)
            if chrom not in coverage:
                coverage[chrom] = {}
            for i in range(int(start), int(end)):
                coverage[chrom][i] = coverage[chrom].get(i,float(0)) + ( count*scale )


if coverage:
    write_bedgraph_from_dict(coverage, args.OUTPUT, digits=2, parallel=False, multi_key=False)



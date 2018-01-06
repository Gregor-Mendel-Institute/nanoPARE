import argparse
import time
import multiprocessing as mp

desc = (
    "Takes a BEDGRAPH of count values.\n"
    "All contiguous regions with values surpassing a given threshold\n"
    "are converted to features in an output BED file.\n"
    "\n"
)

parser = argparse.ArgumentParser(description=desc)

# add arguments to the ArgumentParser
parser.add_argument(
    '-B', '--bedgraph', dest='bedgraph_in', metavar='"file"',
    type=str, help='input bedgraph', action='store',
    required=True
)
parser.add_argument(
    '-O', '--output', dest='bed_out', metavar='"file"',
    type=str, help='output bedgraph', action='store',
    required=True
)
parser.add_argument(
    '-L', '--lengths', dest='lengths', metavar='"file"',
    type=str, help='lengths table', action='store',
    required=True
)
parser.add_argument(
    '-T', '--threshold', dest='threshold', type=float,
    help='minimum threshold above which to call a feature',
    default=0
)
parser.add_argument(
    '-M', '--minimum', dest='minimum', type=int,
    help='minimum length for a feature.',
    default=10
)
parser.add_argument(
    '-V', '--value', dest='value', type=str,
    help='type of score to output in bed file',
    choices=['sum','mean','max'], default='sum'
)
parser.add_argument(
    '-S', '--strand', dest='strand', type=str,
    help='Strand of features to output in bed file',
    choices=['+','plus','p','.','-','minus','m'], default='.'
)

args = parser.parse_args()

if args.strand in ['plus','p']:
    args.strand = '+'

if args.strand in ['minus','m']:
    args.strand = '-'

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]

def find_peaks(vector):
    """Returns local maxima in an array"""
    peaks       = []
    dec         = False
    inc         = False
    plateau     = False
    was_inc     = False
    was_plateau = False
    prev        = 0
    for i in range(len(vector)):
        curr = vector[i]
        if curr > prev:
            inc = True
            dec = False
        elif curr < prev:
            inc = False
            dec = True
        elif curr == prev:
            inc = False
            dec = False
        
        if dec and ( was_inc or plateau ):
            peaks += [i-1]
        
        plateau     = not dec and not inc and ( was_inc or was_plateau )
        was_plateau = plateau
        was_inc     = inc
        prev        = curr
    return peaks

def writer(merge_filename,shared_queue,stop_token):
    """Initializes a writer to output the mp queue."""
    dest_file = open(merge_filename,'w')
    while True:
        line = shared_queue.get()
        if line == stop_token:
            dest_file.close()
            return
        dest_file.write(line)

def find_bed_features(chrom,chromlen,queue):
    above_thresh = sorted(list(thresh_coverage[chrom]))
    if len(above_thresh) < 2:
        return
    
    # Identify discontinuous regions in the positions above the threshold
    splits = [(a,b) for a,b in zip(above_thresh[:-1],above_thresh[1:]) if b - a > 1]
    
    # Pair start and end locations for each region
    region_endpoints = zip([above_thresh[0]] + [b for a,b in splits],[a for a,b in splits] + [above_thresh[-1]])
    
    featurecount = 0
    # Iterate through regions, writing a BED entry if large enough
    for chromStart,chromEnd in region_endpoints:
        if chromEnd - chromStart >= args.minimum:
            featurecount += 1
            name = '.'.join(['thresh',chrom,str(featurecount)])
            coverage_subset = [coverage[chrom].get(i,0) for i in range(chromStart,chromEnd+1)]
            if args.value == 'sum':
                signalValue = sum(coverage_subset)
            elif args.value == 'mean':
                signalValue = sum(coverage_subset)/len(coverage_subset)
            elif args.value == 'max':
                signalValue = max(coverage_subset)
            
            peak_positions = find_peaks(coverage_subset)
            if len(peak_positions) == 0:
                print('WARNING: no peaks found in {}:{}-{}; skipping.'.format(chrom,chromStart,chromEnd))
                continue
            peak_heights = [coverage_subset[i] for i in peak_positions]
            ordered_peaks = [
                peak_positions[i]
                for i in [
                    k
                    for v,k in sorted(
                        [
                            (v,k)
                            for k,v in enumerate(peak_heights)
                        ], reverse=True
                    )
                ]
            ]
            peak = ordered_peaks[0]
            if len(ordered_peaks) > 1:
                other_peaks = [str(i) for i in ordered_peaks[1:]]
            else:
                other_peaks = ['-']
            
            queue.put(
                '\t'.join(
                    [
                        str(i)
                        for i in [
                            chrom,
                            chromStart,
                            chromEnd,
                            name,
                            round(signalValue,2),
                            args.strand,
                            peak,
                            ','.join(other_peaks)
                        ]
                    ]
                ) + '\n'
            )

# 'chromosomes' contains the lengths of all chromosomes 
# that BEDGRAPH contains values for.
# Expects a two-column tab-separated file with:
#     chromosome  length
# Provided with the -L argument.
chromosomes={}
lengths_file=open(args.lengths)
for line in lengths_file:
    chrom,length=line.rstrip().split('\t')
    chromosomes[chrom]=int(length)

# 'coverage' is a dictionary of float vectors 
# for each nucleotide in the genome.
# Contains the value of the BEDGRAPH file at each position.
coverage={}
thresh_coverage={}
for chrom,chromlen in chromosomes.items():
    coverage[chrom] = {}
    thresh_coverage[chrom] = set()

coverage_file = open(args.bedgraph_in)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    for i in range(int(start),int(end)):
        coverage[chrom][i] = count
        if count > args.threshold:
            thresh_coverage[chrom].add(i)

# Begins identifying continuous features, outputting a BED file.
# Example line:
#    Ath_chr1 1242364 1243615 thresh.358 0 + 31.1916 -1 -1 118
#        chrom           chromosome name
#        chromStart      leftmost end of the read (0-indexed)
#        chromEnd        rightmost end of the read (0-indexed, open)
#        name            unique name of each peak
#        score           signal in peak, measured as value type
#        strand          directionality of peak (+ or -)
#        peak            position of peak point source (relative, 0-indexed)

if __name__ == '__main__':
    queue = mp.Queue()
    STOP_TOKEN = "FINISHED"
    writer_process = mp.Process(
        target=writer,
        args=(args.bed_out,queue,STOP_TOKEN)
    )
    writer_process.start()
    all_threads = []
    chromosomes_by_size = [
        k for v,k in sorted(
            [(v,k) for k,v in chromosomes.items()]
            ,reverse=True
        )
    ]
    for chrom in chromosomes_by_size:
            all_threads.append(
                mp.Process(
                    target=find_bed_features,
                    args=(
                        chrom,
                        chromosomes[chrom],
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
    

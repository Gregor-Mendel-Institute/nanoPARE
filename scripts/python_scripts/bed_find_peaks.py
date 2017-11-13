import argparse
import time
import multiprocessing as mp

desc = (
    "Takes a BED file of features\n"
    "and one or more bedgraph files.\n"
    "The highest-coverage position in each feature is output,\n"
    "as well as secondary peaks, relative to the start position.\n"
)

parser = argparse.ArgumentParser(description=desc)

# add arguments to the ArgumentParser
parser.add_argument(
    '-BP', '--bedgraph_plus', dest='bedgraph_plus',
    type=str, help='input bedgraph', nargs='+',
    default=[]
)
parser.add_argument(
    '-BM', '--bedgraph_minus', dest='bedgraph_minus',
    type=str, help='input bedgraph', nargs='+',
    default=[]
)
parser.add_argument(
    '-I', '--input', dest='bed_in', metavar='"file"',
    type=str, help='input bed', action='store',
    required=True
)
parser.add_argument(
    '-O', '--output', dest='bed_out', metavar='"file"',
    type=str, help='output bed', action='store',
    required=True
)
parser.add_argument(
    '-L', '--lengths', dest='lengths', metavar='"file"',
    type=str, help='lengths table', action='store',
    required=True
)

parser.add_argument(
    '-V', '--value', dest='value', type=str,
    help='type of score to output in bed file',
    choices=['pass','sum','rpm','mean','max','winsor','period','byrank','allvalues'], default='pass'
)

args = parser.parse_args()

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

def fold_vector_from_peak(vector):
    """Returns an array folded in half from the peak position
    Position 0 is the highest value, and each position n
    is the sum of the elements +-n from the peak.
    """
    # peak position in vector
    mp = min(which(vector,max(vector)))
    folded_length = max([mp,len(vector)-mp+1])
    folded_vector = [0]*folded_length
    folded_vector[0] = vector[mp]
    for i in range(1,folded_length):
        folded_value = 0
        if mp + i < len(vector):
            folded_value += vector[mp + i]
        if mp - i > 0:
            folded_value += vector[mp - i]
        folded_vector[i] = folded_value
        
    # cut off any zero values from the far end of the vector
    folded_vector = folded_vector[:max(which([i>0 for i in folded_vector]))+1]
    return folded_vector

    
def find_max_period(vector, winsorize=False, quantize=False):
    """Determine which phase contains the highest mean read counts"""
    if winsorize:
        # set a cutoff of 90% of the total score
        target = sum(vector)*.9
        sv = sorted(vector,reverse=True)
        cumsum = [sum(sv[:(i+1)]) for i in range(len(sv))]
        # Find the minimum number of positions to explain >=90% of vector
        score = min(which([float(i) >= target for i in cumsum]))+1
        # Filter out all values less than the height threshold to 0
        filterval = sv[score-1]
        filtervector = [i if i >= filterval else 0 for i in vector]
        vector = filtervector
    else:
        score = len(which([i > 0 for i in vector]))
    
    if quantize:
        vector = [1 if i > 0 else 0 for i in vector]
    
    # Report the most enriched phase in distance from the peak
    # Format: number,phase,proportion in phase
    s = sum(vector)
    if s == 0:
        score = '0,0,0'
    else:
        fv = fold_vector_from_peak(vector)
        if len(fv) <= 1:
            score = '1,1,1'
        else:
            pdict = {}
            for p in range(1,len(fv)):
                in_frame = which([i%p==0 for i in range(len(fv))])
                reads_in_frame = float(sum([fv[i] for i in in_frame]))
                
                meanval = reads_in_frame/len(in_frame)
                proportion_in_frame = reads_in_frame/s
                
                pdict[p] = (meanval*proportion_in_frame,proportion_in_frame)
            
            sorted_results = sorted([(v,k) for k,v in pdict.items()],reverse=True)
            if len(sorted_results) >= 1:
                max_result = sorted_results[0]
                max_p = max_result[1]
                max_meanval = round(max_result[0][0],3)
                max_prop = round(max_result[0][1],3)
            else:
                print('WARNING: {} yielded no results'.format('.'.join([chrom,str(counter)])))
                max_p=max_meanval=max_prop=0
            
            
            score = "{},{},{}".format(score,max_p,max_prop)
    
    return score


def sorted_positions(vector):
    """Returns local maxima in an array"""
    peaks = [k for v,k in sorted([(v,k) for k,v in enumerate(vector)],reverse=True)]
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

def bed_find_peaks(chrom, queue, secondary=False, continuous=False):
    counter=1
    if chrom not in bed_dict:
        return
    
    for i in bed_dict[chrom]:
        start,end,strand,name,score,other = i
        vector = coverage[strand][chrom][start:end]
        if continuous:
            peaks = find_peaks(vector)
        else:
            peaks = sorted_positions(vector)
        if len(peaks) == 0:
            dompeak = '.'
            otherpeaks = '.'
        else:
            dompeak = peaks[0]
            if len(peaks) > 1 and secondary:
                otherpeaks = ','.join([str(j) for j in peaks[1:]])
            else:
                otherpeaks = '.'
        if args.value == 'pass':
            score = score
        elif args.value == 'winsor':
            # set a cutoff of 90% of the total score
            target = sum(vector)*.9
            sv = sorted(vector,reverse=True)
            cumsum = [sum(sv[:(i+1)]) for i in range(len(sv))]
            # Find the minimum number of positions to explain >=90% of vector
            score = min(which([float(i) >= target for i in cumsum]))+1
            # Filter out all values less than the height threshold to 0
            filterval = sv[score-1]
            filtervector = [i if i >= filterval else 0 for i in vector]            
            # Calculate the mean distance between positions with values
            if score <= 1:
                meansep = 0
            else:
                seps = 0
                lastpos = None
                positions = 0
                for pos in range(len(filtervector)):
                    if filtervector[pos]:
                        if lastpos:
                            seps += pos - lastpos
                        lastpos = pos
                        positions += 1
                
                meansep = round(float(seps) / (positions - 1), 3)
            
            score = "{},{}".format(score,meansep)
        elif args.value == 'period':
            score = find_max_period(vector)
        elif args.value == 'byrank':
            # set a cutoff of 90% of the total score
            target = sum(vector)*.9
            sv = sorted(vector,reverse=True)
            cumsum = [sum(sv[:(i+1)]) for i in range(len(sv))]
            # Find the minimum number of positions to explain >=90% of vector
            score = min(which([float(i) >= target for i in cumsum]))+1
            # Filter out all values less than the height threshold to 0
            filterval = sv[score-1]
            filtervector = [i if i >= filterval else 0 for i in vector]            
            # Sort the filtered peak locations by abundance and return their positions as a list
            poslist = [k for v,k in sorted([(v,k) for k,v in enumerate(filtervector)],reverse=True) if v > 0]
            score = ','.join([str(i) for i in poslist])
        elif args.value == 'allvalues':
            s = sum(vector)
            threshold = 0.01
            vector_norm = [float(i)/s for i in vector]
            vector_norm = [round(i,2) if i >= threshold else 0 for i in vector_norm]
            scorevector = []
            for i in range(len(vector_norm)):
                if vector_norm[i]:
                    scorevector += ['{}:{}'.format(i,vector_norm[i])]
            
            if len(scorevector) == 0:
                score = '0:0'
            else:
                score = ','.join(scorevector)
        elif args.value == 'rpm':
            score = round(float(sum(vector))/readnumber*1000000,2)
        else:
            score = sum(vector)
        
        queue.put(
            '\t'.join(
                [
                    chrom,
                    str(start),
                    str(end),
                    '.'.join([chrom,str(counter)]),
                    str(score),
                    strand,
                    str(dompeak),
                    other
                ]
            ) + '\n'
        )
        counter += 1


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
coverage['+'] = {}
coverage['-'] = {}
for chrom,chromlen in chromosomes.items():
    coverage['+'][chrom] = [0]*chromlen
    coverage['-'][chrom] = [0]*chromlen

for file in args.bedgraph_plus:
    coverage_file = open(file)
    for line in coverage_file:
        chrom,start,end,count = line.rstrip().split()
        count = float(count)
        for i in range(int(start),int(end)):
            coverage['+'][chrom][i] += count
    coverage_file.close()

for file in args.bedgraph_minus:
    coverage_file = open(file)
    for line in coverage_file:
        chrom,start,end,count = line.rstrip().split()
        count = float(count)
        for i in range(int(start),int(end)):
            coverage['-'][chrom][i] += count
    coverage_file.close()

readnumber = int(
    sum([sum(i) for i in coverage['+'].values()]) +
    sum([sum(i) for i in coverage['+'].values()])
)

# Imports input bed file, organizing by chromosome
bed_dict = {}
file = open(args.bed_in)
for line in file:
    l = line.rstrip().split('\t')
    chrom = l[0]
    start = int(l[1])
    end = int(l[2])
    name = l[3]
    score = l[4]
    strand = l[5]
    other = l[6]
    if chrom not in bed_dict:
        bed_dict[chrom] = []
    bed_dict[chrom] += [(start,end,strand,name,score,other)]
    
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
                    target=bed_find_peaks,
                    args=(
                        chrom,
                        queue,
                        False
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
    

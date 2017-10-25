import re
import sys
import argparse
import gff_utils as gu
import fasta_utils as fu
import multiprocessing as mp

########################
### ARGUMENT PARSING ###
########################

desc = (
    "Takes a genome-mapped BEDGRAPH file and an annotation file (GFF),"
    "outputs a new BEDGRAPH file where 'chromosomes' are each transcript."
)

parser = argparse.ArgumentParser()
parser.add_argument("bedgraph_plus",
                    help="BEDGRAPH of plus-stranded nucleotide-level coverage.")
parser.add_argument("bedgraph_minus",
                    help="BEDGRAPH of minus-stranded nucleotide-level coverage.")
parser.add_argument(dest='reference_GFF',
                    help="Path to a reference annotation GTF/GFF.")
parser.add_argument(dest='genome_fasta',
                    help="Path to a FASTA file of the genome.")
parser.add_argument('--subset', default='', type=str,
                    help="Path to file of transcript IDs to use from reference.")
parser.add_argument('-O','--output', default='transcript_coverage.bedgraph', type=str,
                    help="Filename to output the transcript-level bedgraph.")
parser.add_argument('--write_fasta', default='', type=str,
                    help="(optional) Filename to output the transcript FASTA file.")


args = parser.parse_args()

####################
# DEFINE FUNCTIONS #
####################

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]


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


####################
# LOAD ENVIRONMENT #
####################

# TODO: retool gff_utils to import gff3 format files
ref_transcripts = {}
refGFF = open(args.reference_GFF)
for line in refGFF:
    if line[0] == '#':continue
    chrom,source,gtype,start,end,score,strand,phase,other = line.rstrip().split('\t')
    if gtype not in ['exon','pseudogenic_exon','miRNA_primary_transcript']:
        continue
    
    parent = re.search('^.*;?Parent=([^;]+);?.*$',other).groups()[0]
    ID = re.search('id=(.+?)[\.:;].+$',other)
    if ID:
        ID = ID.groups()[0]
    else:
        ID = parent.split('\.')[0]
    
    iso_length = len(re.findall(',',parent))+1
    isos = re.search(','.join([ID+'(\.[0-9]+)?']*iso_length),parent).groups()
    iso_ID = parent.split(',')
    for i in iso_ID:
        ref_transcripts[i] = ref_transcripts.get(i,{})
        if gtype != 'miRNA_primary_transcript':
            exon_number = re.search(
                'ID=' + ID + ':(pseudogenic_)?exon:([0-9]+).*;Parent.+$',
                other
            )
            if exon_number:
                exon_number = exon_number.groups()[1]
            else:
                exon_number = str(max([int(j) for j in ref_transcripts[i].get('exon_nums',['0'])])+1)
        else:
            exon_number = '1'
    
    for i in isos:
        if i:
            iso_ID = ID + i
        else:
            iso_ID = ID
        
        ref_transcripts[iso_ID]['chrom'] = chrom
        ref_transcripts[iso_ID]['strand'] = strand
        ref_transcripts[iso_ID]['start'] = ref_transcripts[iso_ID].get('start',[])+[int(start)]
        ref_transcripts[iso_ID]['end'] = ref_transcripts[iso_ID].get('end',[])+[int(end)]
        ref_transcripts[iso_ID]['exon_nums'] = ref_transcripts[iso_ID].get('exon_nums',[])+[exon_number]  
        
ref_IDs = sorted(list(ref_transcripts.keys()))
print('# {} reference transcripts: {}'.format(
    len(ref_IDs),
    args.reference_GFF
))

# 'picked_IDs' is an array of IDs to use from the reference_GFF
if args.subset:
    picked_IDs = [i.rstrip() for i in open(args.subset).readlines()]
else:
    picked_IDs = ref_IDs

# 'genome' is a dict of strings for each chromosome in 'genome_fasta'
genome = fu.import_genome(args.genome_fasta)

# 'chromosomes' contains the lengths of all chromosomes the that BEDGRAPH contains values for.
chromosomes = {}
for chrom in genome.keys():
    length = len(genome[chrom])
    chromosomes[chrom] = int(length)

# 'coverage' is a dictionary of float vectors for each nucleotide in the genome.
# Contains the value of the BEDGRAPH file at each position.
coverage = {}
coverage['plus'] = {}
coverage['minus'] = {}

for chrom,chromlen in chromosomes.items():
    coverage['plus'][chrom] = [0]*chromlen
    coverage['minus'][chrom] = [0]*chromlen

# Populate the plus bedgraph dict
coverage_file = open(args.bedgraph_plus)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    coverage['plus'][chrom][int(start):int(end)] = \
        [count]*(int(end)-int(start))

# Populate the minus bedgraph dict
coverage_file = open(args.bedgraph_minus)
for line in coverage_file:
    chrom,start,end,count = line.rstrip().split()
    count = float(count)
    coverage['minus'][chrom][int(start):int(end)] = \
        [count]*(int(end)-int(start))

# Open a FASTA file for the picked IDs
if args.write_fasta:
    fasta_outfile = open(args.write_fasta,'w')

# Generate an output dict of transcript-level coverages
output_coverage_dict = {}

for ID in picked_IDs:
    if ID not in ref_transcripts:
        print('# WARNING: {} not found'.format(ID))
        continue
    
    chrom = ref_transcripts[ID]['chrom']
    strand = ref_transcripts[ID]['strand']
    exon_starts = sorted(ref_transcripts[ID]['start'])
    exon_ends = sorted(ref_transcripts[ID]['end'])
    
    positions = flatten(
        [list(range(a - 1,b)) for a,b in zip(exon_starts,exon_ends)]
    )
    
    if strand == '+':
        transcript_coverage = [coverage['plus'][chrom][i] for i in positions]
        if args.write_fasta:
            transcript_fasta = ''.join([genome[chrom][i] for i in positions])
            fasta_outfile.write('>{}\n'.format(ID))
            fasta_outfile.write('{}\n'.format(transcript_fasta))
    elif strand == '-':
        transcript_coverage = list(reversed([coverage['minus'][chrom][i] for i in positions]))
        if args.write_fasta:
            transcript_fasta = fu.rc(''.join([genome[chrom][i] for i in positions]))
            fasta_outfile.write('>{}\n'.format(ID))
            fasta_outfile.write('{}\n'.format(transcript_fasta))
    else:
        print('# WARNING: {} is unstranded'.format(ID))
        continue
    
    output_coverage_dict[ID] = {}
    for i in range(len(transcript_coverage)):
        output_coverage_dict[ID][i] = transcript_coverage[i]

if args.write_fasta:
    fasta_outfile.close()

# Write transcript-level bedgraph to the output file
write_bedgraph_from_dict(
    output_coverage_dict,
    output_filename=args.output
)


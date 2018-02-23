# Used to estimate the optimal scaling ratio between a
# sequencing library of end reads (5'/3') against a
# sequencing library of body reads (RNA-seq)

import re
import sys
import argparse

######################################################################
parser = argparse.ArgumentParser()
parser.add_argument('END_QUANT',
    help="Path to a positive quant table (assumes reads/length/RPM/TPM)."
)
parser.add_argument('BODY_QUANT',
    help="Path to a background quant table (assumes reads/length/RPM/TPM)."
)
parser.add_argument('--baseline', dest='BASELINE',
    help="Manually set a baseline ratio for end/body fragment number.",
    default=None, type=float
)
args = parser.parse_args()
######################################################################

def parse_quant_table(quant_table):
    file = open(quant_table)
    quant = {}
    for line in file:
        if line[0] == '#':
            continue
        
        l = line.rstrip().split('\t')
        gene,reads,length,rpm,tpm = l
        quant[gene] = (float(reads),int(length),float(rpm),float(tpm))
    return quant

if args.BASELINE:
    baseline = args.BASELINE
else:
    #TODO: calculate an estimated baseline empirically
    baseline = 0.25

end_quant = parse_quant_table(args.END_QUANT)
body_quant = parse_quant_table(args.BODY_QUANT)

# Get a list of genes shared in common between the two tables
end_genes = set(list(end_quant.keys()))
body_genes = set(list(body_quant.keys()))
shared_genes = sorted(list(end_genes.intersection(body_genes)))

# Count the total number of reads mapping to all shared genes
reads_end = 0
reads_body = 0
for g in shared_genes:
    reads_end += end_quant[g][0]
    reads_body += body_quant[g][0]

# Print a scaling factor representing the value of 1 end read
# relative to 1 body read, given the expected baseline
print(round(reads_body/(reads_end/baseline),3))

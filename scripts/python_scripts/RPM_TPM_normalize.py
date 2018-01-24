import argparse

########################
### ARGUMENT PARSING ###
########################

parser = argparse.ArgumentParser()
parser.add_argument('-R','--rpm',dest='RPM',
                    help="Path to a RPM table.",
                    required=True)
parser.add_argument('-T','--tpm',dest='TPM',
                    help="Path to TPM table.",
                    required=True)
parser.add_argument('-L','--length',dest='FRAGMENT_LENGTH',
                    help="Mean fragment length of the sample.",
                    default=200)
args = parser.parse_args()

# Calculate a conversion ratio to normalize the sequencing depth of
# and end-specific library with that of a full-length library

# total_effective_length is the sum total length (in nucleotides)
# of a representative collection of 1 million transcripts laid end-to-end
total_effective_length = 0
total_BODY_reads = float(0)
total_END_reads = float(0)

for line in open(args.TPM):
    if line[0] == '#':
        continue
    l = line.rstrip().split()
    raw = float(l[1])
    length = int(l[2])
    norm = float(l[3])
    if raw > 0 and norm > 0:
        total_BODY_reads += raw
        total_effective_length += length*norm

# BODY_population is the total number of fragments
# that would be generated from a pool of 1 million transcripts
BODY_population = total_effective_length/10**6/args.FRAGMENT_LENGTH
END_BODY_ratio = 1 / BODY_population

for line in open(args.RPM):
    if line[0] == '#':
        continue
    l = line.rstrip().split()
    raw = float(l[1])
    length = int(l[2])
    norm = float(l[3])
    if raw > 0 and norm > 0:
        total_END_reads += raw

print( END_BODY_ratio / (total_END_reads / total_BODY_reads))

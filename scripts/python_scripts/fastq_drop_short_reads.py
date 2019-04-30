import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--minlen", dest='MINLEN',
    help="Minimum length for a FASTQ read to keep.",
    default=18, type=int
)
parser.add_argument(
    "-O", "--output", dest='OUTPUT', type=str, default='stdout',
    help="Filepath to write processed FASTQ file."
)
parser.add_argument(
    "FILENAME", nargs='?'
)
args = parser.parse_args()

####################

if args.FILENAME:
    if args.FILENAME.split('.')[-1].lower() not in ['fq','fastq']:
        print("\nERROR: input file must be FASTQ format.")
        parser.print_help()
        sys.exit(1)
    fastq_in = open(args.FILENAME)
elif not sys.stdin.isatty():
    fastq_in = sys.stdin
else:
    print("\nERROR: requires FASTQ file as input.")
    parser.print_help()
    sys.exit(1)

if args.OUTPUT != 'stdout':
    args.OUTPUT = open(args.OUTPUT,'w')

####################

def output_bed_lines(output_lines,output=args.OUTPUT):
    """Takes a list of bed lines and writes
    them to the output stream.
    """
    if output == 'stdout':
        for output_string in output_lines:
            print(output_string)
    else:
        for output_string in output_lines:
            output.write('{}\n'.format(output_string))

entry1=[]
linecounter=0
tooshortcount=0
for i in fastq_in:
    linecounter+=1
    line1=i.rstrip()
    if linecounter % 4 == 1:
        entry1=[]
    entry1.append(line1)
    if linecounter % 4 == 0:
        if len(entry1[1])>= args.MINLEN:
            output_bed_lines(entry1,args.OUTPUT)
        else:
            tooshortcount+=1

if args.OUTPUT != 'stdout':
    args.OUTPUT.close()

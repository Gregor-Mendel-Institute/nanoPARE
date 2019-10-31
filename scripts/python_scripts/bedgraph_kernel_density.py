#!/usr/bin/env python3
__author__ == "Falko Hofmann"
__email__ == "falkohofmann@gmail.com"
"""
TODO: add module docstring
"""
import logging
import argparse
import time
from argparse import RawTextHelpFormatter
from collections import OrderedDict
import multiprocessing as mp
from multiprocessing import Pool
import numpy as np

LOGGER = logging.getLogger('bedgraph_kernel_density')
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s')

def parse_args():
    """Function that returns an argparse object for parsing the required
    parameters.

    Returns:
        parser: Reuturns parser.parse_args() object

    Raises:
        argparse.ArgumentTypeError
    """
    # strings explaining the usage, displayed by argparse
    desc = ("Takes a BEDGRAPH of count values."
            " Each point is converted to a density distribution and all"
            " densities are summed for the output model.")

    epilog = ("Any values immediately upstream of -U or downstream of -D will"
              " be set to zero.This can be used to mask sequence-specific"
              " artifacts, like TSO strand invasion or oligo-dT mispriming.")

    # initilize argumentparser to read in commands
    parser = argparse.ArgumentParser(description=desc,
                                     epilog=epilog,
                                     formatter_class=RawTextHelpFormatter)

    # add arguments to the ArgumentParser
    parser.add_argument('-B', '--bedgraph', dest='bedgraph_in', action='store',
                        type=str, help='input bedgraph', metavar='"file"',
                        required=True)
    parser.add_argument('-O', '--output', dest='bedgraph_out', action='store',
                        type=str, help='output bedgraph', metavar='"file"',
                        required=True)
    parser.add_argument('-L', '--lengths', dest='lengths', metavar='"file"',
                        type=str, help='lengths table', action='store',
                        required=True)
    parser.add_argument('-K', '--kernel', dest='kernel', metavar='"str"',
                        help='Type of the kernel function.', type=str,
                        default='laplace', choices=['gaussian', 'laplace'])
    parser.add_argument('-H', '--bandwidth', dest='bandwidth', metavar='int',
                        type=int, help='Bandwith of the kernel', default=15)
    parser.add_argument('-S', '--sigma', dest='sigma', metavar='integer',
                        type=int, help='sigma value', choices=range(1, 7),
                        default=3)
    parser.add_argument('-D', '--digits', dest='digits', metavar='int',
                        type=check_positive, help='Name of the output file',
                        default=3)
    parser.add_argument('-P', '--positive', dest='positive', default=False,
                        action='store_true', help='Write only positive values')
    parser.add_argument('--single_stream', dest='single_stream', default=False,
                        action='store_true', help='(for very large files) process bedgraph data as a single stream')
    parser.add_argument('-c', '--cores', dest='cores', metavar='int',
                        type=int, help='Number of CPU cores to use.',
                        default=1)

    return parser.parse_args()  # parse everything


# function for argparse arugment checking
def check_positive(value):
    """Function that checks if a variable passed on to argparse is positive
       integer'.

    Args:
        value: Variable to check

    Returns:
        int: The input variable as int.

    Raises:
        argparse.ArgumentTypeError
    """
    # TODO write pytest with strings as input
    i = int(value)
    if i < 0:
        raise argparse.ArgumentTypeError('%s is an invalid positive int value',
                                         value)
    return i


def mp_writer_dict(dest_dict,shared_queue,stop_token):
    """Initializes a writer to output the mp queue."""
    while True:
        line = shared_queue.get()
        if line == stop_token:
            return
        tkey, tval = line
        dest_dict[tkey] = tval


def read_file(f):
    """Function that returns all lines contained in a the file 'f'.

    Args:
        f (str): A file path

    Returns:
        list: The lines of the file.
    """
    logging.info('Parsing file %s', f)
    with open(f) as fin:
        return fin.readlines()


def get_laplace_kernel(bandwidth, sigma):
    """Function that returns a simplified laplace kernel

    Args:
        bandwidth(int): Bandwidth of the kernel
        sigma(int): Sigma of the kernel

    Returns:
        k(numpy.array): An array representation of the kernel
    """
    logging.info('Generating laplace kernel with bandwidth: %s, sigma: %s',
                 bandwidth, sigma)
    r = np.arange(-bandwidth*sigma, bandwidth*sigma+1, dtype=np.float32)
    k = laplace(r, bandwidth)
    return normalize_kernel(k)


def get_gaussian_kernel(bandwidth, sigma):
    """Function that returns a simplified gaussian kernel

    Args:
        bandwidth(int): Bandwidth of the kernel
        sigma(int): Sigma of the kernel

    Returns:
        k(numpy.array): An array representation of the kernel
    """
    logging.info('Generating gaussian kernel with bandwidth: %s, sigma %s',
                 bandwidth, sigma)
    r = np.arange(-bandwidth*sigma, bandwidth*sigma+1, dtype=np.float32)
    k = gauss(r, bandwidth)
    return normalize_kernel(k)


def normalize_kernel(kernel):
    """Function that scales a kernel to a sum of 1.

    Args:
        kernel(numpy.array): The kernel to scale.

    Returns:
        k(numpy.array): A numpy array scaled to a sum of 1.
    """
    return kernel/np.sum(kernel)


def gauss(x, s):
    """Function that calulates a gaussian distribution for a numpy array'.
       The calculation performed is: e**(-x**2/2*s**2)/sqrt(2*pi)*s

    Args:
        x(numpy array): Values across which the gaussian distribution should be
                        calculated
        s(int): Sigma value for the gaussian distribution.

    Returns:
        numpy array: The transformed values.
    """
    # calculate the first exponent of e
    e_exp_x = -1*np.power(x, 2)
    # calculate the second exponent of e
    e_exp_s = np.power(2*s, 2)
    # divide the exponents
    e_exp = np.divide(e_exp_x, e_exp_s)
    # calculate the nominator
    nom = np.exp(e_exp)
    # then the denominator
    denom = np.sqrt(2*np.pi) * s
    # divide and return
    return np.divide(nom, denom)


def laplace(x, s):
    """Function that calulates a simplified laplace distribution (m = 0),
       for a numpy array'. The calculation performed is: e**(-(|x|-m)/s)/2*s

    Args:
        x(numpy array): Values across which the laplace distribution should be
                        calculated
        s(int): Sigma value for the laplace distribution.
    Returns:
        A numpy array of the transformed values.
    """
    # TODO maybe add implementation of muy when this will be supplied as
    # stand-alone package
    # x_neg = np.absolute(x) * -1
    # calculate nominator
    nom = np.exp(np.divide(np.absolute(x), s)*-1)
    # then the denominator
    denom = np.array([2*s])
    # and divide and return
    return np.divide(nom, denom)
    # return math.exp(-abs(float(x)-m)/s)/(2*s)


def parse_chromosomes(file_name):
    """
    Function that parses a file that contains the length of all chromosomes.
    Provided with the -L argument.

    Args:
        file_name(str): A two-column tab-separated file with: chromosome length

    Returns:
        A dict mapping chromosome names to their length. For example:
        {'Chr1': 30427671, 'Chr2': 19698289}
    """
    # Use an ordered dictionary to allow ordered output writing.
    d = OrderedDict()
    lines = read_file(file_name)
    for l in lines:
        # FIXME skipp emtpy lines
        chrom, length = l.rstrip().split('\t')
        length = int(length)
        # some sanity check, chromomes need to be larger than >1
        if length < 1:
            logging.warning('Chromosome %s has invalid length with %s.'
                            'Skipping.', chrom, length)

            continue
        else:
            logging.info('Chromosome %s has length  %s.', chrom, length)
            d[chrom] = length
    return d


def init_chromosomes(chroms):
    """
    Function that initializes a dictionary of numpy arrays for each nucleotide
    in the genome.

    Args:
        d_chrom(dict): A dictionary mapping chromosome names ot their length.
        As e.g. the output from parse_chromosomes().

    Returns:
        A dict mapping chromosome names to a numpy array of their length with
        zero as default value.
    """
    # Use an ordered dictionary to allow ordered output writing.
    d_cov = OrderedDict()
    for chrom, chromlen in chroms.items():
        logging.info('Initializing chromosome %s with %s.', chrom, chromlen)
        d_cov[chrom] = np.zeros(chromlen)
    
    return d_cov


def get_bedgraph_coverage(chroms, file_name):
    """
    Function that reads in a bedgraphs file and returns the values.

    Args:
        chroms(dict): A dictionary mapping chromosome names ot their length.
        file_name(str): Path of the file to parse.

    Returns:
        A dict mapping chromosome names to a numpy array. In the numpy array
        each position reflects a base in the genome and the bedgraph value.
    """
    cov = init_chromosomes(chroms)
    lines = read_file(file_name)
    for l in lines:
        chrom, start, stop, value = l.rstrip().split('\t')
        # dont try to put data into non initialized chromosomes
        if chrom not in cov:
            logging.warning('Length of chromosome %s has not been specified. '
                            'Skipping bedgraph region.', chrom)
            continue
        start = int(start)
        stop = int(stop)
        value = float(value)
        # do some sanity checks on the input data to avoid complications
        too_big = any(x > len(cov[chrom]) for x in (start, stop))
        too_small = any(x < 0 for x in (start, stop))
        if too_big or too_small:
            logging.warning('Bedgraph contains malformed entry. Skipping.\n'
                            '%s', l)
            continue
        # put data after succcessfull sanity checks
        cov[chrom][start:stop] = value
    
    return cov


def get_pool_args(d, k, pos):
    # FIXME: Method not needed anymore
    """
    Function that pares a list of parameters for do_convolution().

    Args:
        d(dict): A dict mapping chromosome_id to numpy.array
        k(numpy.array): Kernel that should later be convoluted on d.values()
        pos(bool): Boolean that indicates of only positive values should be
                   returned after kernel smoothing.

    Returns:
        Returns a list of tuples of the structure
        [(chrom_id, coverage, kernel, return_only_above_zero)]
    """
    l = []
    # As memory is not shared accross processes, they all need their instance
    # of data. This loops helping putting it into the right format.
    for key, val in d.items():
        # tuples are easier to feed into methods run with Pool
        l.append((key, val, k, pos))
    return l


def do_convolution(data, parallel=False, queue=None):
    """
    Function that performs a convoution of kernel k on coverage data from a
    specific chromosome. This function is designed to run with multiprocessing
    Pool().

    Args:
        data(tuple): A tuple containing ((chromosome_id(str),
        coverage(numpy.array), kernel(numpy.array), above_zero(bool))

    Returns:
        Returns a tuple of the structure (chromosome_id(str),
        smoothed_coverage (numpy.array)). Depending on the boolean flag
        'above_zero' the array contains only values above 0 or all values.
    """
    chrom, cov, k, only_positive = data
    k_midpoint = int((len(k)-1)/2)
    position_offsets = [i - k_midpoint for i in range(len(k))]
    # logging.info('Fitting kernel to chromosome %s:', chrom)
    print ('Fitting kernel to chromosome %s.' % chrom)
    # Do a convolution with a kernel along the coverage to smoothen it
    cov_smoothed = np.convolve(cov, k, 'same')
    if only_positive:
        neg_values = np.where(cov_smoothed < 0)
        cov_smoothed[neg_values] = 0
    
    if parallel:
        queue.put([chrom,cov_smoothed])
    else:
        return chrom, cov_smoothed


def get_changes(cov, same_values, compare_zero):
    """
    Generator function that checks if and how the value in an array is changing
    beween the i and i-1 position the array.

    Args:
        arr(list): A list contining a coverage profile

    Yields:
        Returns a tuple containing the type of the change (e.g. from_zero)
        the position and the value. Only yields non zero regions.
    """
    start = 0
    stop = 0
    for i, v in enumerate(cov):
        stop = i + 1
        if compare_zero[i]:
            start = i
            continue
        if not same_values[i]:
            temp = start
            start = i + 1
            yield(temp, stop, v)


def coverage_to_bedgraph(data):
    """
    Function that converts an array containing coverage data into a bedgraph
    like format.

    Args:
        data(tuple): A tuple of 3 contining chromosome id, coverage and number
        of significant digits.

    Returns:
        Returns a list of tuples(id, start, stop, value) containing the
        bedgraph bins.
    """
    name, cov, digits = data
    cov = np.round(cov, digits)
    # logging.info('Converting chromosome %s to bedgraph.', name)
    print('Converting chromosome %s to bedgraph.' % name)
    is_zero = np.isclose(cov, 0)
    cov_c = np.roll(cov, -1)
    cov_c[-1] = 0
    is_same = np.isclose(cov, cov_c)
    bedgraph = []
    # begraph is a region/bin format, so convert the continuos arrays to bins
    # with non 0 values
    for start, stop, value in get_changes(cov, is_same, is_zero):
        bedgraph.append((name, start, stop, value))
    return bedgraph


def single_stream_process(input_file,output_file,chromosomes,kernel,digits=2,only_positive=True):
    """
    Perform all operations for kernel density estimation on a 
    single bedgraph file as it is being read. Minimizes RAM requirements
    for processing very large files
    """
    def is_not_zero(value,only_positive=False):
        """ Returns whether a value is greater than zero, non-zero, or zero """
        if only_positive:
            return value > 0
        else:
            return value != 0
    
    
    def dump_chunk(chunk,chrom,startpos,output,only_positive):
        """ Writes all elements of the chunk """
        start = startpos
        prevpos = 0
        prevcount = None
        position = None
        
        for current_value in chunk:
            if position is None:
                position = start
            else:
                position += 1
            
            if position < 0:
                continue
            
            count = round(current_value,digits)
            if count != prevcount or int(position) > 1 + prevpos:
                # The newly encountered value is not a continuation
                # of the previous value. Write the old run and start another.
                if prevcount and is_not_zero(prevcount,only_positive):
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
                    output.write(line_to_write)
                        
                start = position
            prevcount = count
            prevpos = position
        
        if position and prevcount and is_not_zero(prevcount,only_positive):
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
            output.write(line_to_write)
    
    
    input = open(input_file)
    output = open(output_file,'w')
    # Get kernel size and relative positions
    klen = len(kernel)
    k_midpoint = int((klen-1)/2)
    position_offset = -k_midpoint
    current_chrom = None
    # Initialize a dict that will temporarily store positions within the zone of effect of the kernel
    current_chunk = []
    chunk_start = 0
    for l in input:
        # Perform all calculations on a local chunk
        chrom, start, stop, value = l.rstrip().split('\t')
        # dont try to put data into non initialized chromosomes
        if not current_chrom:
            current_chrom = chrom
            print ('Fitting kernel to chromosome %s.' % chrom)
        
        if chrom != current_chrom:
            print ('Fitting kernel to chromosome %s.' % chrom)
            dump_chunk(current_chunk, current_chrom, chunk_start, output, only_positive)
            current_chunk = []
            chunk_start = 0
            current_chrom = chrom
        
        start = int(start)
        stop = int(stop)
        value = float(value)
        if value == 0:
            continue
        
        # Perform kernel convolution with the current line
        for pos in range(start,stop):
            kval = [i*value for i in kernel]
            kpos = pos + position_offset
            if kpos < chunk_start:
                # Expand the left edge of chunk to accommodate the kernel
                current_chunk = [0]*(chunk_start - kpos) + current_chunk
                chunk_start = kpos
            kernel_endpos = pos - position_offset
            chunk_endpos = chunk_start + len(current_chunk)
            if kernel_endpos > chunk_endpos:
                if kpos > chunk_endpos:
                    # Dump the chunk and start a new one
                    dump_chunk(current_chunk, current_chrom, chunk_start, output, only_positive)
                    current_chunk = [0]*klen
                    chunk_start = kpos
                    chunk_endpos = chunk_start + klen
                else:
                    # Expand the right edge of chunk to accommodate the kernel
                    current_chunk = current_chunk + [0]*(kernel_endpos - chunk_endpos)
            # Add the kernel to current_chunk
            for pos_to_add,val_to_add in zip(range(kpos - chunk_start,kernel_endpos - chunk_start),kval):
                current_chunk[pos_to_add] += val_to_add
    
    input.close()
    dump_chunk(current_chunk, current_chrom, chunk_start, output, only_positive)
    output.close()



def main():
    """
    TODO: Add docstring.
    """
    # parse all the arguments to determine what should be done
    args = parse_args()
    bedgraph_in = args.bedgraph_in
    bedgraph_out = args.bedgraph_out
    lengths = args.lengths
    kernel_type = args.kernel
    bandwidth = args.bandwidth
    sigma = args.sigma
    digits = args.digits
    only_positive = args.positive
    nprocs = args.cores
    single_stream = args.single_stream
    parallel = nprocs > 1
    # calculate the kernel so that it can be convolved over the data
    kernel = None
    if kernel_type == 'laplace':
        kernel = get_laplace_kernel(bandwidth, sigma)
    elif kernel_type == 'gaussian':
        kernel = get_gaussian_kernel(bandwidth, sigma)
    # read in the length of the chromosomes
    chromosomes = parse_chromosomes(lengths)
    if single_stream:
        # Process bedgraph data as single input/output streams
        single_stream_process(
            bedgraph_in,
            bedgraph_out,
            chromosomes,
            kernel,
            digits,
            only_positive
        )
    else:
        # read in the begraph data so that it later can be smoothed
        coverage = get_bedgraph_coverage(chromosomes, bedgraph_in)
        # restructure some objects in a list of tuples that can be used easily
        # with multiprocessing
        convolution_args = get_pool_args(coverage, kernel, only_positive)
        # Do the convolution with the selected kernel to get a smoothed coverage
        results = []
        for data in convolution_args:
            results.append(do_convolution(data))
        # reshape the smoothed coverage into a tuple that can be used to construct
        # a bedgraph
        # IDEA: switch from digit rounding before comparison to only printing the
        # respective numbers
        data_out = [r + (digits,) for r in results]
        # convert the smoothed coverage data to a bedgraph like format
        # via multiprocees and write to output file.
        p = Pool(nprocs)
        with open(bedgraph_out, 'w') as f:
            for results in p.imap(coverage_to_bedgraph, data_out):
                for name, start, stop, value in results:
                    fmt_str = '%s\t%d\t%d\t%%.%sf\n' % (name, start, stop, digits)
                    f.write(fmt_str % value)


# run everything
if __name__ == '__main__':
    main()

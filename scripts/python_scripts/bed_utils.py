#!/usr/bin/env python3
__author__ == "Falko Hofmann"
__email__ == "falkohofmann@gmail.com"

"""
This  module is for parsing of bed files.
"""
# TODO: add parsing safety checks.

import logging
from genomic_features import Transcript

LOGGER = logging.getLogger('bed_utils')
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s')

(CHROM_IDX, START_IDX, STOP_IDX, ID_IDX, SCORE_IDX, STRAND_IDX, THICKSTART_IDX,
 THICKEND_IDX, RGB_IDX, BLOCK_COUNT_IDX, BLOCK_SIZE_IDX, BLOCK_START_IDX) = range(0, 12)

def read_files(paths):
    """Method that reads a bed file and returns its content

    Args:
        paths (list: str): A list of bed file paths that should be parsed.

    Returns:
        regions (dict): A dictionary mapping the region name to a tuple of
            (region name(str), start(int), stop(int)).

    """
    LOGGER.info('Reading %d files...', len(paths))
    content = {}
    for path in paths:
        content.update(read_file(path))
    LOGGER.info('Reading %d files... - Done', len(paths))
    return content


def read_file(path):
    """Method that reads a bed file and returns its content

    Args:
        path (str): Path of a bed file that should be parsed.

    Returns:
        regions (dict): A dictionary mapping the region name to a tuple of
            (region name(str), start(int), stop(int)).

    """
    regions = {}
    with open(path) as fin:
        LOGGER.info('Reading bed file: %s', path)
        for name, region in process_lines(fin):
            regions[name] = region
        LOGGER.info('Reading bed file: %s - Done', path)
    return regions


def process_lines(lines):
    """Generator function that splits the columns of a bed file.

    Args:
        lines (list): A list of bed file lines.

    Yields:
        results (tuple: str, int, int): A tuple of the bed region id, start and
            stop position

    """

    for line in lines:
        cols = line.split('\t')
        results = (cols[ID_IDX], (int(cols[START_IDX]), int(cols[STOP_IDX])))
        yield results

def read_bed12(path, source, sample_name, mode='dictionary'):
    transcripts = {}
    with open(path) as fin:
        LOGGER.info('Reading bed file: %s', path)
        for name, transcript in process_bed12_lines(fin, source, sample_name):
            transcripts[name] = transcript
        LOGGER.info('Reading bed file: %s - Done', path)
    return transcripts


def process_bed12_lines(lines, source, sample_name):
    for line in lines:
        bed_cols = line.split('\t')
        transcript = Transcript(*bed_cols, formatting='bed12', source=source)
        transcript.add_sample_name(sample_name)
        yield transcript.name, transcript

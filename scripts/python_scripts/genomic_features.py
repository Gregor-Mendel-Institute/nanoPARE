#!/usr/bin/env python3
"""
This module is a data container for transcripts parsed via gff_utils or bed_utils.
Transcripts should be the only class you actually need to import.

Example:
    from genomic_features import Transcripts

"""
__author__ = "Falko Hofmann"
__email__ = "falkohofmann@gmail.com"

import logging
from abc import ABC, abstractmethod
from collections import namedtuple

#supported and allowed file formats
ALLOWED_FORMATS = frozenset(['bed12', 'gtf/gff'])
# file type column headers
BED12_COLS = ['chromosome', 'start', 'stop', 'name', 'score', 'strand',
              'thickstart', 'thickstop', 'rgb', 'block_count', 'block_size',
              'block_start']
GTFGFF_COLS = ['chromosome', 'source', 'feature', 'start', 'stop', 'score',
               'strand', 'frame', 'attributes']
#create some named tuples to abstract out the index numbers
Bed12Values = namedtuple('Bed12Values', BED12_COLS)
BED12 = Bed12Values(*range(0, len(BED12_COLS)))
GtfgffValues = namedtuple('GtfgffValues', GTFGFF_COLS)
GTFGFF = GtfgffValues(*range(0, len(GTFGFF_COLS)))

#create logger
LOGGER = logging.getLogger('genomic_features')
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s')

class Feature(ABC):
    """A class abstraction of a gff3/gtf genomic feature. Inherits from ABC.
    Implements most of the attributes of a genomic feature that are different
    for each type of feature. Class Transript and Exon inherit from it.

    Args:
        data (list): Columns of gff3/gtf line
        feature_id (str): A string holding a unique identifier.

    """

    # TODO: move some of the attributes to class specific attritbues
    def __init__(self, feature_id, formatting, *args, attributes=None):
        # initilize variables and fill with some data
        # TODO: recheck which variables are needed
        self.name = None
        self.score = None
        self.strand = None
        self.phase = None
        self.start = None
        self.stop = None
        self.attributes = None

        self.name = feature_id

        if formatting == 'bed12':
            self.score = args[BED12.score]
            self.strand = args[BED12.strand]
            self.phase = '.'
            self.start = args[BED12.start]
            self.stop = args[BED12.stop]
            self.attributes = attributes

        if formatting == 'gtf/gff':
            self.score = args[GTFGFF.score]
            self.strand = args[GTFGFF.strand]
            self.phase = args[GTFGFF.frame]
            self.start = args[GTFGFF.start]
            self.stop = args[GTFGFF.stop]
            self.attributes = args[GTFGFF.attributes]

        self.length = self.stop - self.start

    def update_attributes(self, attributes):
        """Setter method for updating transcript attributes.

        Args:
            attributes (dict): Dictionary that contains key, value pairs that
                should be added as attributes.

        """
        self.attributes = {**self.attributes, **attributes}

    def get_attributes(self, file_type):
        """Method that generates a string from a key,value dictionary containing
        the genomic feature attributes (Column 9) in either as gtf or gff3 format.
        Can be used for converting between file formats and printing.

        Args:
            file_type (str): String indicating format that should be returned.
                Needs to be either 'gff3' or 'gtf'.

        Returns:
            attributes (str): Feature attributes (Column 9 in gff3/gtf) as a string
                in gff3 or gtf format.

        """
        # The attribute colum consists out of key value pairs sitting in
        # dictionary. This loop glues them back together with the appropriate
        # file delimiter

        # TODO: add raise exception.
        attribs = []
        for key, value in self.attributes.items():
            if file_type == 'gtf':
                attribs.append(' '.join((key, '"' + value + '"')))
            if file_type == 'gff3':
                attribs.append('='.join((key, value)))
        if file_type == 'gtf':
            return '; '.join(attribs) + ';'
        if file_type == 'gff3':
            return ';'.join(attribs)

    def set_start_stop(self, start, stop):
        """Setter method for start and stop attributes. Requires start < stop.

        Args:
            start (int): The start position.
            stop (int): The stop position.

        Raises:
            AttributeError: Will be raised when start > stop.

        """
        if start > stop:
            raise AttributeError('Start position of s% is larger than it\'s'
                                 ' stop position.', self.name)
        self.start = start
        self.stop = stop

    @abstractmethod
    def set_feature_id(self, feature_id):
        """ Setter method to set the feature id. Classes inherting from
        feature can implemented their own method.

        Args:
            feature_id (str): A unique id used to describe this transcript.

        """
        self.name = feature_id
        self.attributes['transcript_id'] = feature_id

    @abstractmethod
    def get_feature_id(self, attributes):
        """Getter method to return the unique feature id. Needs to be
        implemented by classes inherting from Feature.

        """
        pass

    @abstractmethod
    def sanitize(self):
        """Method to clear attributes, phase and score of the feature object.
        """
        self.attributes = {}
        self.phase = '.'
        self.score = '.'


class Transcript(Feature):
    """A class abstraction of a gff3/gtf transcript . Inherits from Feature.
    Acts as container for data and Exon objects.

    Args:
        data (list): Columns of gff3/gtf line containing transcript
            information

    """
    # IDEA: add attribute for merged or not
    # TODO: add attribtue for references transcript or not
    # TODO: make parser read alternate tss/tes sites and their signal
    # TODO: clean up constructor so it works with either GTF or BED12
    def __init__(self, *args, formatting='gtf/gff', source=None):

        if formatting not in ALLOWED_FORMATS:
            raise AttributeError('Format must be either %s', ALLOWED_FORMATS)

        self.chrom = None
        self.source = None
        self.exons = []
        self.exon_order = []
        self.number_exons = 0
        self.sample_names = set()
        self.alt_start = []
        self.alt_stop = []
        self.tss_signal = []
        self.tes_signal = []
        self.summed_signal = []
        self.ref_transcript = False

        args = list(args)
        if args and formatting == 'gtf/gff':
            args[GTFGFF.start] = args[GTFGFF.start]
            args[GTFGFF.stop] = args[GTFGFF.stop]
            super().__init__(self.get_feature_id(args[GTFGFF.attributes]),
                             formatting, *args)
            self.chrom = args[GTFGFF.chromosome]
            self.source = args[GTFGFF.source]

        elif args and formatting == 'bed12':
            # bed is 0 based index, gtf/gff 1
            # IDEA maybe move casting of ints to another class for consistency
            args[BED12.start] = int(args[BED12.start]) + 1
            # args[BED12.stop] = int(args[BED12.stop]) + 1
            args[BED12.stop] = int(args[BED12.stop])
            attributes = self.__make_bed12_attributes(args[BED12.name])
            super().__init__(args[BED12.name], formatting, *args, attributes=attributes)
            self.chrom = args[BED12.chromosome]
            self.source = source

            for data in self.__make_exon_ranges(args):
                self.add_exon(data, formatting=formatting, attributes=attributes)

    def __repr__(self):
        #id_str_1 = ''.join(map(str, self.get_exon_start()))
        #id_str_2 = ''.join(map(str, self.get_exon_end()))
        #return self.strand + self.chrom + id_str_1 + id_str_2 + self.name
        #return self.get_bio_id() + ':' + str(self.start) + '_' + str(self.stop)
        return self.name

    def __make_exon_ranges(self, args):
        data = args.copy()
        exons_start = args[BED12.block_start]
        exon_length = args[BED12.block_size]
        exon_start = map(int, exons_start.split(','))
        exon_length = map(int, exon_length.split(','))
        for start_pos, length in zip(exon_start, exon_length):
            exon_start = self.start + start_pos
            exon_stop = exon_start + length - 1
            data = args.copy()
            data[BED12.start] = exon_start
            data[BED12.stop] = exon_stop
            yield data

    def __make_bed12_attributes(self, name):
        transcript_id = '{name}.{iso}'.format(name=name, iso=1)
        return {'gene_id': name, 'transcript_id': transcript_id}

    def get_difference(self, other):
        diff_start = max(self.start, other.start) - min(self.start, other.start)
        diff_stop = max(self.stop, other.stop) - min(self.stop, other.stop)
        diff_total = diff_start + diff_stop

        diff_exons = None
        if self.number_exons == 1 and other.number_exons == 1:
            diff_exons = 0
        elif self.number_exons == 1:
            if self.can_be_merged_with(other, 'any'):
                diff_exons = other.number_exons - self.number_exons
            else:
                diff_exons = other.number_exons
        elif other.number_exons == 1:
            if other.can_be_merged_with(self, 'any'):
                diff_exons = self.number_exons - other.number_exons
            else:
                diff_exons = self.number_exons
        else:
            exons_self = self.get_exon_start()
            exons_other = other.get_exon_start()
            set_self = set()
            set_self.add(exons_self[1])
            set_self.add(exons_self[-2])
            for idx in range(3, self.number_exons - 1, 2):
                start = exons_self[idx]
                stop = exons_self[idx + 1]
                set_self.add(start + stop)

            set_other = set()
            set_other.add(exons_self[1])
            set_other.add(exons_self[-2])
            for idx in range(3, other.number_exons - 1, 2):
                start = exons_other[idx]
                stop = exons_other[idx + 1]
                set_other.add(start + stop)

            common = set_self & set_other
            diff_exons = max(self.number_exons, other.number_exons) - len(common)

        return diff_exons, diff_total


    def set_reference_transcript(self, value):
        self.ref_transcript = value

    def is_reference_transcript(self):
        return self.ref_transcript

    def set_start_sites(self, start_sites, storted_sites=None):
        """Setter method for alternative transcript start sites. Overwrites
        existing values.

        Args:
            start_sites (list: int): List of the alternative start sites.
            sorted_sites (list: int): Arbitrarily sorted version of start_sites.
                This list will be converted to a string (exluding the primary
                start site of the transcript), which is used for neater
                display in the gtf/gff output.

        """
        self.alt_start = start_sites
        # TODO: maybe add some attibute setter when sorted sites are not provided
        if storted_sites:
            if self.strand == '+':
                attibute_key = 'alternative_5p'
            if self.strand == '-':
                attibute_key = 'alternative_3p'
            self.attributes[attibute_key] = self.uniquify_to_string(storted_sites,
                                                                    self.start)

    def set_stop_sites(self, stop_sites, storted_sites=None):
        """Setter method for alternative transcript stop sites. Overwrites
        existing values.

        Args:
            stop_sites (list: int): List of the alternative stop sites.
            sorted_sites (list: int): Arbitrarily sorted version of stop_sites.
                This list will be converted to a string (exluding the primary
                stop site of the transcript), which is used for neater
                display in the gtf/gff output.
        """
        self.alt_stop = stop_sites
        # TODO: maybe add some attibute setter when sorted sites are not provided
        if storted_sites:
            if self.strand == '+':
                attibute_key = 'alternative_3p'
            if self.strand == '-':
                attibute_key = 'alternative_5p'
            self.attributes[attibute_key] = self.uniquify_to_string(storted_sites,
                                                                    self.stop)

    def uniquify_to_string(self, to_uniquify, to_filter):
        """Method to generate a string of unique values from an iterable.

        Args:
            to_uniquify (iterable): A collection of elments that can be iterated over.
            to_filter (Object): A value that should be filtred out of 'to_uniquify'.
                The object needs to be hashable and comparable.
        Returns:
            unique_seq (str): A string of the input elements, with duplicates and
                to_filter removed. Elements are seperated by ','. Returns 'NA' if
                to_uniquify is None.
        """
        # if not to_uniquify:
        #     return "NA"
        seen = set([to_filter])
        seen_add = seen.add
        unique_seq = [str(x) for x in to_uniquify if not (x in seen or seen_add(x))]
        if not unique_seq:
            return "NA"
        return ','.join(unique_seq)



    def get_sample_names(self):
        """Getter method for the sample names.

        Returns:
            sample_names (set): A set of the sample names this transcript is
                associatied with.

        """
        return self.sample_names

    def get_start_sites(self):
        """Getter method for transcript start sites.

        Returns:
            start_sites (list: int): List containing all known transcript start
                sites.

        """
        return self.alt_start if self.alt_start else [self.start]

    def get_stop_sites(self):
        """Getter method for transcript stop sites.

        Returns:
            stop_sites (list: int): List containing all known transcript stop
                sites.

        """
        return self.alt_stop if self.alt_stop else [self.stop]


    def get_lowest_start_site(self):
        """Getter method that returns the lowest transcript start site.

        Returns:
            start_site (int): Lowest start site set for this transcript.

        """
        return min(self.alt_start) if self.alt_start else self.start

    def get_highest_start_site(self):
        """Getter method that returns the highest transcript start site.

        Returns:
            start_site (int): Highest start site set for this transcript.

        """
        return max(self.alt_start) if self.alt_start else self.start

    def get_lowest_stop_site(self):
        """Getter method that returns the lowest transcript stop site.

        Returns:
            stop_site (int): Lowest stop site set for this transcript.

        """
        return min(self.alt_stop) if self.alt_stop else self.stop

    def get_highest_stop_site(self):
        """Getter method that returns the highest transcript stop site.

        Returns:
            stop_site (int): Highest stop site set for this transcript.

        """
        return max(self.alt_stop) if self.alt_stop else self.stop

    def set_tss_signal(self, signal):
        """Setter method for the tss_signal attribute. Overwrites existing
        values.

        Args:
            signal (list: float): A list of floats.

        """
        self.tss_signal = signal

    def set_tes_signal(self, signal):
        """Setter method for the tes_signal attribute. Overwrites existing
        values.

        Args:
            signal (list: float): A list of floats.

        """

        self.tes_signal = signal


    def get_tss_signal(self):
        """Getter method for the signal at the transcript start site (TSS).

        Returns:
            signal (float): The signal as float.

        """

        if 'TSS' in self.attributes:
            return [float(self.attributes['TSS'])]
        if not self.tss_signal:
            return[0]
        return self.tss_signal
        #    raise KeyError('No TSS signal set for this transcript.')

    def get_tes_signal(self):
        """Getter method for the signal at the transcript end site (TES).

        Returns:
            signal (float): The signal as float.

        """

        if 'TES' in self.attributes:
            return [float(self.attributes['TES'])]
        if not self.tes_signal:
            return [0]
        return self.tes_signal
        #else:
        #    raise KeyError('No TSS signal set for this transcript.')

    def get_total_signal(self):
        """Getter method for the signal at the transcript end site (TES).

        Returns:
            signal (float): The signal as float.

        """

        return [sum(x) for x in zip(self.get_tss_signal(),
                                    self.get_tes_signal())]

    def get_tss_id(self):
        """Getter method fot the transcription start site peak id.

        Returns:
            tss_id (str): The id of the tss peak that is associatied with this
                transcript.

        Raises:
            KeyError: A key error will be raised if the tss id has not been set.

        """

        if 'TSS_ID' in self.attributes:
            return self.attributes['TSS_ID']
        raise KeyError('No TSS ID set for this transcript.')

    def get_tes_id(self):
        """Getter method fot the transcription end site peak id.

        Returns:
            te`s_id (str): The id of the tes peak that is associatied with this
                transcript.

        Raises:
            KeyError: A key error will be raised if the tes id has not been set.

        """

        if 'TES_ID' in self.attributes:
            return self.attributes['TES_ID']
        raise KeyError('No TSS ID set for this transcript.')

    def add_sample_name(self, sample_name):
        """Setter method that adds a sample_name to the transcript attributes.

        Args:
            sample_name (str): The sample name that should be added.

        """
        self.sample_names.add(sample_name)

    def update_sample_names(self, sample_names):
        """Setter method that adds multiple sample_names to the transcript
        attributes.

        Args:
            sample_names (list: str): The sample names that should be added.

        """

        self.sample_names.update(sample_names)
        self.attributes['built_in'] = ','.join(sorted(self.sample_names))
        # TODO: Maybe add warning when stuff gets overwritten.

    def set_start_stop(self, start, stop):
        """Setter method for start and stop attributes. Requires start < stop.

        Args:
            start (int): The start position.
            stop (int): The stop position.

        Raises:
            AttributeError: Will be raised when start > stop.

        """
        # no need to for a special function implementation, just call
        # the one from the super module and make sure the exons are properly
        # updated.
        super().set_start_stop(start, stop)
        self.update_exons()
        #self.update_exon_boundaries()

    def update_exons(self):
        """Method that checks if the exons contained in a transcript are
        compatible with the transcript start and stop positions and updates
        the class attributes accordingly. When there are exons that lie outside
        of the transcript start stop boundaries they will be ommitted.
        """

        temp_exons = self.exons
        self.exons = []
        for exon in temp_exons:
            # outside of new transcript
            if exon.start < self.start and exon.stop < self.start:
                continue
            # outside of new transcript
            if exon.start > self.stop and exon.stop > self.stop:
                continue
            if exon.start <= self.start < exon.stop:
                exon.set_start(self.start)
            if exon.start < self.stop <= exon.stop:
                exon.set_stop(self.stop)
            self.exons.append(exon)
        # FIXME: fix exon numbers
        self.exon_order = self.sort_exons()
        self.number_exons = len(self.exons)
        if not self.exons:
            LOGGER.warning('Updated transcript contains no exons %s', self.name)
    # def update_exon_boundaries(self):
    #     # if self.strand == '+'
    #     self.exons[self.exon_order[0]].set_start(self.start)
    #     self.exons[self.exon_order[-1]].set_stop(self.stop)
    #
    #     for x in [0, -1]: # Loop over the first and last exon
    #         if self.exons[x].start > self.exons[x].stop:
    #             LOGGER.warning('Start position of %s is larger than it\'s stop '
    #                            'position. Start: %d Stop: %d', self.name,
    #                            self.exons[x].start, self.exons[x].stop)
    # elif self.strand == '-':
    #     self.exons[self.exon_order[-1]].set_start(self.start)
    #     self.exons[self.exon_order[0]].set_stop(self.stop)
    # else:
    #     raise AttributeError('Strand for %s is set to %s. Cannot update exon boundaries when the strand is '
    #                          'neither + or -.', self.name, self.strand )
    def set_feature_id(self, feature_id):
        """Class specific setter implementation of the method. Sets the
        feature id also for all the exons.

        Args:
            feature_id (str): A unique id used to describe this transcript.

        """

        super().set_feature_id(feature_id)
        for exon in self.exons:
            exon.set_feature_id(feature_id)

    def get_bio_id(self):
        """Method that returns an id that is a summary of the strand,
        chromosome and inner splice junction positions. This id can be used to
        group similar transcripts for biological meaningful comparisons (e.g.
        transcripts that share the same strand, chromosome and splice
        junctions). However comparing Transcript objects should be done
        by the '=='/'!=' operators or the compare_intron_chain method as this
        is slightly faster.

        Returns:
            id (str): An id string summarizing the transcript properties.
                Consists out of strand + chromosome + exon_start_positions
                + exon_stop_positions without the outer end positions (e.g.
                +Ath_chr198457989089780598605).

        """

        if self.number_exons != 1:
            id_str_1 = '_'.join(map(str, self.get_exon_start()[1:]))
            id_str_2 = '_'.join(map(str, self.get_exon_end()[:-1]))
            return self.strand + self.chrom + '_' + id_str_1 + '_' + id_str_2
        # FIXME not an ideal solution. Need to come up with a way how to
        # deal with single exon transcripts
        return (self.strand + self.chrom + str(self.start) + str(self.stop) +
                '_1_exon')

    # self = (29741736, 29743598)
    # other = (29741772, 29743190)
    # --> false should be true
    # self = (29741772, 29743190)
    # other = (29741736, 29743598)
    # --> true


    def get_exon_hash(self):
        # TODO: check if this function is still needed.
        """Method that returns a hash that is a summary of the intron exon
        junctions.

        Returns (str): Hash comprised out exon1Start_exon1Stop_exon2Start etc.

        """

        exon_starts = map(str, self.get_exon_start()[1:])
        exon_stops = map(str, self.get_exon_end()[:-1])
        exon_hash = []
        for positions in zip(exon_starts, exon_stops):
            exon_hash.append('_'.join(positions))
        return '_'.join(exon_hash)

    def overlaps_with(self, other):
        """Method that checks if one transcript overlaps with another, or if
        one is contained in another one.

        Args:
            other (Transcript): The other transcript to compare against.

        Returns:
            check (boolean): Pass of fail of the check.

        """

        # overlap is not possible if the transcripts are on other strands or chromosomes
        if self.strand != other.strand:
            return False
        if self.chrom != other.chrom:
            return False

        # overlaps can either happen on the left or right boundary
        left_overlap = (other.start <= self.start <= other.stop)
        right_overlap = (other.start <= self.stop <= other.stop)

        # is the one contained in the other? which is also an overlap..
        is_contained = (self.start >= other.start and
                        self.stop <= other.stop)
        contains = (other.start >= self.start and
                    other.stop <= self.stop)

        return left_overlap or right_overlap or is_contained or contains

    def compare_intron_chain(self, transcript):
        """Method to compare the intron chain of this transcript with that from
        another. This method only compares if the junctions share the same
        positions. It does not require exact matches of 5'and 3' ends of the
        fist and last exon. Returns true when both transcripts have no introns.

        Args:
            transcript (Transcript): An object of the class Transcript.

        Returns:
            check (boolean): Boolean indicating a failed or succcessfull check.

        """

        if self.chrom != transcript.chrom:
            return False
        elif self.strand != transcript.strand:
            return False

        exons_self = set(self.get_exon_start()[1:])
        exons_self.update(self.get_exon_end()[:-1])

        exons_other = set(transcript.get_exon_start()[1:])
        exons_other.update(transcript.get_exon_end()[:-1])
        common = exons_self & exons_other

        n_junctions = max(len(exons_self), len(exons_other))
        return len(common) == n_junctions

    def compatible_truncations(self, other):
        """Method to check if one of the transcripts is a truncation of the
        other with a compatible intron chain architecture.

        Args:
            other (Transcript): The other transcript to compare against.

        Returns:
            compatible (boolean): True/False if  one of the transcript a
                truncation with compatible intron chain architecture.

        """

        exons_self = set(self.get_exon_start()[1:])
        exons_self.update(self.get_exon_end()[:-1])
        exons_other = set(other.get_exon_start()[1:])
        exons_other.update(other.get_exon_end()[:-1])

        transcript_start = None
        transcript_stop = None
        intron_ranges = None
        exon_ranges = None

        # exon/intron positions of one need to be subset of the other for a truncation.
        # also do things bidirectional.
        if exons_other.issubset(exons_self):
            intron_ranges = list(zip(self.get_exon_end()[:-1],
                                     self.get_exon_start()[1:]))
            exon_ranges = list(zip(other.get_exon_start(),
                                   other.get_exon_end()))
            # transcript_start = other.get_lowest_start_site()
            # transcript_stop = other.get_highest_stop_site()
            transcript_start = other.start
            transcript_stop = other.stop
        elif exons_self.issubset(exons_other):
            intron_ranges = list(zip(other.get_exon_end()[:-1],
                                     other.get_exon_start()[1:]))
            exon_ranges = list(zip(self.get_exon_start(), self.get_exon_end()))
            #transcript_start = self.get_lowest_start_site()
            #transcript_stop = self.get_highest_stop_site()
            transcript_start = self.start
            transcript_stop = self.stop
        else:
            return False

        # checks if one transcript starts in the intron of another.
        # If this happens, its not a valid truncation.
        for start, stop in intron_ranges:
            if ((start < transcript_start < stop) or
                    (start < transcript_stop < stop)):
                return False
        # checks if intron retention at the terminal exons happens.
        # Not a valid truncation if that happens.
        for e_start, e_stop in exon_ranges:
            for i_start, i_stop in intron_ranges:
                if (e_start < i_start < e_stop) or (e_start < i_stop < e_stop):
                    return False

        return True


    def is_contained_in_exon(self, other):
        """Method that checks if the transcript is fully contained in an exon
        of another transcript. This check is directional and tests if the
        transcript on which this method is called is contained in the one
        provided as function argument.

       Args:
            other (Transcript): The other transcript to compare against.

        Returns:
            compatible (boolean): Returns True/False if the transcript this
                method is called upon is contained in an exon of the 'other'
                transcript.
        """

        #loop over the exon borders of the other transcript one by one and
        # check if it contains this transcript
        ranges = zip(other.get_exon_start(), other.get_exon_end())
        for start, stop in ranges:
            if (start <= self.start <= stop) and (start <= self.stop <= stop):
                return True
        return False

    def terminal_exon_overlap_with(self, other):
        """Checks for terminal exons of one transcript overlapping with another
        single exon transcript. The check requires the other single exon
        transcript to adhere to the intron/exon boundaries of the other. This
        check is directional.

       Args:
            other (Transcript): The transcript to check for terminal overlap with.

        Returns:
            overlap (boolean): Is the 'other' transcript is a single exon
                transcript overlapping in a concordant way with the terminal
                exons?

        """

        # comparison only works and makes sense when one transcript is single exon.
        # the comparison is also directional.
        if other.number_exons != 1:
            return False

        first_exon = self.exons[0]
        last_exon = self.exons[-1]
        #if other.get_highest_stop_site() <= first_exon.stop:
        if other.stop <= first_exon.stop:
            return True
        #if other.get_lowest_start_site() >= last_exon.start:
        if other.start >= last_exon.start:
            return True
        return False


    def can_be_merged_with(self, other, mode):
        """Method that checks if one transcript can be merged with another
        based on some basic biological/set rules. This comparison partially
        abstracts out the directionality and specific cases of the comparisons.
        Depending on the properties of the two transcripts that are beeing
        compared the following checks are performed/returned. Merge modes can
        either set to 'strict' or 'any'.


        'Any' merging behaviour:
            1. Single exon vs. single exon:
                - Checks if the transcripts are overlapping.
            2. Single exon vs multi exon:
                - Is the single exon transcript overlapping with the terminal
                  exons of the multi exon transcript? And does the overlap
                  respect intron/exon boundaries of the multi exon transcript?
                - Or is the single exon transcript completly contained in one
                  exon of the multi exon transcripts?
            3. Multi exon vs multi exon:
                - If both transcripts have the same number of exons, do they
                  ialso have the same intron chain?
                - If they have different number of exons, is one a truncation of
                  the other with the same but incomplete intron/exon structure?

        'Strict' merging behaviour:
            1. Multi exon transcripts:
                - Only returns true for complete intron chain matches.
            2. Single exon transcripts:
                - Only returns true for overlapping transcipts or if one is
                  contained in the other.
            3. Mulit exon vs single exon transcipts:
                - Always returns false

        Args:
            other (Transcript): The other transcript that should be checked.
            scope (str): Needs to be either 'any' or 'strict'.

        Returns:
            check (boolean): Can the 'other' transcript be merged with the one
                this method was called from?

        Raises:
            AttributeError: Will be raised when scope does not equal 'any' or
                'strict'.
        """

        if mode == 'strict':
            if self.number_exons == 1 and other.number_exons == 1:
                return self.overlaps_with(other)
            return self.compare_intron_chain(other)
        elif mode == 'any':
            if self.number_exons == 1 and other.number_exons == 1:
                return self.overlaps_with(other)
            elif other.number_exons == 1:
                return (other.is_contained_in_exon(self) or
                        self.terminal_exon_overlap_with(other))
            elif self.number_exons == 1:
                return (self.is_contained_in_exon(other) or
                        other.terminal_exon_overlap_with(self))
            elif self.number_exons == other.number_exons:
                return self.compare_intron_chain(other)
            return self.compatible_truncations(other)
        else:
            raise AttributeError('Mode needs to be either \'any\' or \'strict\' not  %s', mode)


    def add_exon(self, data, formatting='gtf/gff', attributes=None):
        """Setter method to add exons that belong to a transript.

        Args:
            data (list): Columns of gff3/gtf line containing exon
                information

        """
        # exon_number is needed in case no exon number is provided by the
        # annotation file
        exon_number = len(self.exons) + 1
        exon = Exon(exon_number, *data, formatting=formatting, attributes=attributes)
        self.exons.append(exon)
        self.exon_order = self.sort_exons()
        self.number_exons = len(self.exons)

    def sort_exons(self):
        """Method to sort the exons according to ascending starting position.
        This method will be called after an exon has be added.

        """
        start_l = []
        for exon in self.exons:
            start_l.append(exon.start)
        return sorted(range(len(start_l)), key=start_l.__getitem__)

    def get_exon_start(self):
        """Getter method that returns exon start positions. Exons are sorted
        in an ascending order by their start position.

        Returns:
            positions (list): Ascendingly sorted positions.

        """
        return [self.exons[i].start for i in self.exon_order]

    def get_exon_end(self):
        """Getter method that returns exon end positions. Exons are sorted in
        an ascending order by their start position.

        Returns:
            positions (list): Ascendingly sorted positions.

        """
        return [self.exons[i].stop for i in self.exon_order]

    def get_feature_id(self, attributes):
        """Getter method that returns a unique feature id contained in the
        attributes ('ID' or 'transcript_id').

        Returns:
            feture_id (str): The string of an unique id

        Raises:
            AttributeError: Will be raised when no 'ID' or 'transcript id'
                field is present in the attributes.

        """
        if 'transcript_id' in attributes:
            return attributes['transcript_id']
        elif 'ID' in attributes:
            return attributes['ID']
        elif 'id' in attributes:
            return attributes['id']
        else:
            raise AttributeError('No \'ID\' or \'transcript_id\' attribute'
                                 ' present in the input file. Check your input'
                                 ' file')

    def get_as_string(self, file_type):
        """Method that generates a multi-line string from a transcript and
        its exons in either gff3 or gtf format.

        Args:
            file_type (str): String indicating format that should be returned.
                Needs to be either 'gff3' or 'gtf'.

        Returns:
            attributes (str): A multi-line string in the specified file format.

        """
        lines = []
        first_cols = ('%s\t%s' % (self.chrom, self.source))

        line = ('%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s' %
                (first_cols, 'transcript', self.start, self.stop, self.score,
                 self.strand, self.phase, self.get_attributes(file_type)))

        lines.append(line)

        for exon in self.exons:
            line = ('%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s' %
                    (first_cols, 'exon', exon.start, exon.stop, exon.score,
                     exon.strand, exon.phase, exon.get_attributes(file_type)))
            lines.append(line)

        return '\n'.join(lines) + '\n'


    def sanitize(self):
        """Method to clear attributes, phase and score of the feature object.

        """
        super().sanitize()
        for exon in self.exons:
            exon.sanitize()


class Exon(Feature):
    """A class abstraction of a gff3/gtf exon . Inherits from Feature.
    Acts as container for data.

    Args:
        data (list): Columns of gff3/gtf line containing transcript
            information
        exon_number (int): A number indicating the number of the exon.

    """
    def __init__(self, exon_number, *args, formatting='gtf/gff', attributes=None):
        self.parent = None
        self.exon_number = exon_number
        if args and formatting == 'gtf/gff':
            self.parent = self.get_parent(args[GTFGFF.attributes])
            super().__init__(self.get_feature_id(args[GTFGFF.attributes]), formatting,
                             *args)
        elif args and formatting == 'bed12':
            self.parent = self.get_parent(attributes)
            super().__init__(self.get_feature_id(attributes), formatting, *args, 
                             attributes=attributes)
        else:
            pass
            # TODO: throw error


    def get_parent(self, attributes):
        """Getter method that returns a the parent id (transcript id) of an
        exon.

        Returns:
            feture_id (str): The id of the parent of the exon.

        Raises:
            AttributeError: Will be raised when no 'Parent' or 'transcript id'
                field is present in the attributes.

        """
        if 'Parent' in attributes:
            return attributes['Parent']
        elif 'parent' in attributes:
            return attributes['parent']
        elif 'transcript_id' in attributes:
            return attributes['transcript_id']
        else:
            raise AttributeError('Could not find \'Parent\' or'
                                 ' \'transcript_id\' attribute for an exon.'
                                 ' Check your input file')

    def get_feature_id(self, attributes):
        """Getter method that returns a unique exon id contained in the
        attributes ('ID' or 'exon_id'). If neither of the fields is present
        the method will create an id based on the exon number.

        Returns:
            feture_id (str): The string of an unique id

        """
        if 'ID' in attributes:
            return attributes['ID']
        elif 'id' in attributes:
            return attributes['id']
        elif 'exon_id' in attributes:
            return attributes['exon_id']
        return 'exon:%s:%s' % (self.parent, self.exon_number)


    def set_start(self, start):
        # TODO: add raise error?
        # TODO: make consistent with the other methods
        """Setter method for the start position of an exon.

        Args:
            start (int): The new start position.

        """
        if self.stop < start:
            LOGGER.warning('Start position of %s is larger than it\'s stop position.'
                           ' Start: %d Stop: %d', self.name, start, self.stop)
        self.start = start

    def set_stop(self, stop):
        # TODO: add raise error?
        # TODO: make consistent with the other methods
        """Setter method for the stop position of an exon.

        Args:
            stop (int): The new start position.

        """
        if stop < self.start:
            LOGGER.warning('Start position of %s is larger than it\'s stop position.'
                           ' Start: %d Stop: %d', self.name, self.start, stop)
        self.stop = stop


    def set_feature_id(self, feature_id):
        """Setter method to set the feature id.

        Args:
            feature_id (str): A unique id used to describe this transcript.

        """
        super().set_feature_id(feature_id)
        self.parent = feature_id

    def sanitize(self):
        """Method to clear attributes, phase and score of the feature object.

        """
        super().sanitize()

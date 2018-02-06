import logging
import re
from collections import defaultdict
import sys
if sys.version_info[0] > 2:
    from abc import ABC, abstractmethod
else:
    from abc import ABCMeta, abstractmethod

# import pandas as pd


LOGGER = logging.getLogger('gff_utils')
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s')
# LOGGER.setLevel(20) setting logger level in the module is not recommended

# TODO: check if this is accessible uppon import as module
(CHROM_IDX, SOURCE_IDX, FEATURE_IDX, START_IDX, STOP_IDX, SCORE_IDX,
 STRAND_IDX, FRAME_IDX, ATTR_IDX) = range(0, 9)


def parse_annotation(path, mode='dictionary'):
    """Wrapper function that parses transcripts and exons from an annotation
    file. The function checks the file extention (\'.gtf\' or \'.gff\') and
    then decides how to proceed.

    Args:
        path (str): A file path.

    Returns:
        transcripts: A dictionary or a pandas data frame.

    """
    if path.lower().endswith('.gtf'):
        return parse_gtf(path, mode)
    elif path.lower().endswith('.gff'):
        return parse_gff3(path, mode)
    else:
        raise ValueError('Invalid file type extention! Make sure the file ends'
                         'with either \'.gtf\' or \'.gff\' (Case insensitive).')


def parse_gff3(path, mode='dictionary'):
    """Wrapper function that parses transcripts and exons from a GFF3 file.
    Supports Ensemble, Gencode and TAIR files.

    Args:
        path (str): A file path

    Returns:
        transcripts: A dictionary or a pandas data frame.

    """
    return get_file_content(path, 'gff3', mode)


def parse_gtf(path, mode='dictionary'):
    """Wrapper function that parses transcripts and exons from a GTF file.
    Supports Ensemble, Gencode and TAIR files.

    Args:
        path (str): A file path.

    Returns:
        transcripts: A dictionary or a pandas data frame.

    """
    return get_file_content(path, 'gtf', mode)


def check_file_type_args(file_type):
    """Error handling function. Raises an expetion if a provided string is
    neither 'gtf' or 'gff3'.

    Args:
        file_type (str): A string describing the file_type.

    Raises:
        ValueError

    """
    # TODO add logging
    file_types = {'gtf', 'gff3'}
    if file_type not in file_types:
        raise ValueError('Invalid file type selected! Select either \'gtf\' or'
                         '\'gff3\'')


def get_grouped_transcripts(paths, file_type):
    """Function that processes and returns all transcripts contained in an
    GFF3 or GTF file as a dictionary and groups them by an id key as returned
    by 'transcript.get_bio_id()', which is a a summary of an identical intron
    chain, strand and chromosome. Single exon transcripts will only be grouped
    by strand and chromosome.

    Args:
        paths (iterable): An iterable containing file paths.
        file_type (str): The type of annotation file ('gtf' or 'gff3').

    Returns:
        transcripts (defaultdict(list)): A dictionary mapping where transcripts
            are grouped based on overlapping identiy.
            E.g. key: transcript.get_bio_id(), value: transcript

    """
    # TODO implement auto file detection
    # TODO implement pandas compability
    # if mode not in {'dictionary', 'pandas'}:
    #     raise ValueError('Invalid file mode selected! Select either'
    #                      '\'dictionary\' or \'pandas\'')
    #
    # LOGGER.info('Parsing file: %s Return type: %s', path, mode)
    # raise an error if someone tries to give undefined parameter values
    check_file_type_args(file_type)
    # get regex patterns for appropriate processing according to file type
    patterns = get_regex_pattern(file_type)
    # some containers to hold correct transcripts and orphaned subfeatures
    # which can arise from an unsorted file
    transcripts = defaultdict(list)
    # get the file content
    for path in paths:
        with open(path) as fin:
            # generator to skipp over not needed data and to make things neater
            for transcript in process_lines(fin, patterns):
                # add transcripts and exon features in propper hierarchy
                transcripts[transcript.get_bio_id()].append(transcript)
    return transcripts



def get_file_content(path, file_type, mode='dictionary'):
    """Function that processes and returns all transcripts contained in an
    GFF3 or GTF file either as a dictionary or a pandas data frame.

    Args:
        path (str): A file path
        file_type (str): The type of annotation file ('gtf' or 'gff3').
        mode (str): The return type of the data.

    Returns:
        transcripts (dictionary): A dictionary mapping transcript id to a
            Transcript class object.

    Raises:
        ValueError: Will be raised when the mode is neither 'dictionary' nor
            'pandas'.

    """
    # IDEA: maybe replace with a function that takes a set and a string
    if mode not in {'dictionary', 'pandas'}:
        raise ValueError('Invalid file mode selected! Select either'
                         '\'dictionary\' or \'pandas\'')

    LOGGER.info('Parsing file: %s Return type: %s', path, mode)
    # raise an error if someone tries to give undefined parameter values
    check_file_type_args(file_type)
    # get regex patterns for appropriate processing according to file type
    patterns = get_regex_pattern(file_type)
    # some containers to hold correct transcripts and orphaned subfeatures
    # which can arise from an unsorted file
    transcripts = {}
    # get the file content
    with open(path) as fin:
        # generator to skipp over not needed data and to make things neater
        for transcript in process_lines(fin, patterns):
            # add transcripts and exon features in propper hierarchy
            transcripts[transcript.name] = transcript
    if mode == 'dictionary':
        return transcripts
    if mode == 'pandas':
        # temp_df = pd.DataFrame.from_dict(transcripts, orient='index')
        return temp_df.rename(columns={0:'transcript'})



def process_lines(lines, patterns):
    """Generator processes and yields transcripts and exond from a GFF3 or GTF
    file.

    Args:
        lines (iterable): Lines in the file
        patterns (tuple): A tuple of the compiled regex patterns returned from
            get_regex_pattern()

    Yields:
        transcript (Transcript): An object containing the transcript and exon
            data.

    """
    col_pattern, attr_pattern, field_pattern, field_delimiter = patterns
    transcript = None
    for i, l in enumerate(lines):
        l = l.strip()
        l_num = i + 1

        # only keep transcripts and exons
        is_exon, is_transcript = check_feature_type(l)
        if not any((is_exon, is_transcript)):
            continue
        # regex based column extraction allows some syntax checking
        cols = col_pattern.findall(l)
        # check if the column extraction lead to the expected number (9)
        if not has_9_cols(cols, l_num):
            continue
        # cast start and stop values to integers to check for invalid lenghts
        cols = cast_to_integers(cols)
        # some sanity checks to prevent adding of erroneous lines
        checks = (positive_length(cols, l_num),
                  valid_attributes(cols, field_pattern, field_delimiter, i))
        if not all(checks):
            transcript = None
            continue
        # replace string from the attribute column with key value pairs dict
        cols = process_attributes(cols, attr_pattern)
        if is_transcript:
            if transcript is not None:
                yield transcript
            transcript = Transcript(cols)
        else:
            if transcript is not None:
                if is_exon_of(cols, transcript):
                    transcript.add_exon(cols)
    yield transcript


def is_exon_of(cols, transcript):
    """Function that checks if the columns provided are an exon of the provided
    transcript.
    Args:
        cols (list): A gtf3 or gff line split into different columns.
        transcript (Transcript): An object of the class transcript.

    Returns:
        check (boolean): If the transcript name is contained within exon
            attributes.

    """
    return transcript.name in cols[ATTR_IDX].values()

def cast_to_integers(cols):
    """Function that casts the start (column 4) and stop (column 5) position of
    a genomic feature from strings into integers.

    Args:
        cols (list): A gtf3 or gff line split into different columns

    Returns:
        cols (list): Updated cols, with start and stop columns cast to integers.

    """
    # IDEA maybe add try catch?
    cols[START_IDX] = int(cols[START_IDX])
    cols[STOP_IDX] = int(cols[STOP_IDX])
    return cols


def positive_length(cols, line_num):
    """Function that checks the length validity of the genomic feature.

    Args:
        cols(list): A gtf3 or gff line split into different columns
        line_number(int): Number of the line to check for loggin purposes.

    Returns:
        check (boolean): Boolean indicating a failed or succcessfull check.

        """
    length = cols[STOP_IDX] - cols[START_IDX]
    check = length > 0
    if not check:
        LOGGER.warning('Line: %d contains invalid start and stop positions',
                       line_num)
    return check


def process_attributes(cols, pattern):
    """Function that processes the gff3/gtf attributes column (column 9) into a
        dictionary mappin {'attribute name': 'value'}

    Args:
        cols (list): A gtf3 or gff line split into different columns
        pattern: Compiled regex pattern returned from get_regex_pattern()
            to extract attribute values.

    Returns:
        cols (list): Updated cols, with the attribute column restructured into a
            dictionary.

    """
    # TODO check possible forbidden special chars and add them to the pattern
    attr = cols[ATTR_IDX]
    i = iter(pattern.findall(attr))
    d_attr = dict(zip(i, i))
    cols[ATTR_IDX] = d_attr
    return cols


def has_9_cols(cols, line_num):
    """Function that checks if a line has the expected number of columns.

    Args:
        cols (list): A gtf3 or gff line split into different columns
        line_number(int): Number of the line to check for loggin purposes.

    Returns:
        check (boolean): Boolean indicating a failed or succcessfull check.

    """
    check = len(cols) == 9
    if not check:
        LOGGER.warning('Line: %d contains an invalid number columns', line_num)
    return check


def check_feature_type(line):
    """Function that checks if a line from a gff3/gtf file is either a
    transcript or an exon.

    Args:
        line (str): An unproceessed gtf or gff3 line.

    Returns:
        tuple (boolean): Tuple indicating if the line is an exon or a transcript

    """
    return is_exon(line), is_transcript(line)


def is_exon(line):
    """Function that checks if a gff3/gtf line represents an exon.

    Args:
        line (str): An unproceessed gtf or gff3 line.

    Returns:
        check (boolean): Boolean indicating if the line contains an exon or not.

    """
    return '\texon\t' in line.lower()


def is_transcript(line):
    """Function that checks if a gff3/gtf line represents a transcript.

    Args:
        line (str): An unproceessed gtf or gff3 line.

    Returns:
        check (boolean): Boolean indicating if the line contains a transcript

    """
    # TODO check if this are all possible words denoting a transcript
    check_set = {'\ttranscript\t', '\tmrna\t'}
    checks = [s in line.lower() for s in check_set]
    return any(checks)


def get_regex_pattern(file_type):
    """Function returns a set of regex patterns to split and process either a
    gff3 or gtf formatted line.

    Args:
        file_type (str): String idicating the file type. Has to be either 'gff3'
            or 'gtf'.

    Returns:
        col_pattern (re.pattern): Compiled regex pattern to split a line into
            the number of appropriate columns. Imposes some mild character
            restrictions. To ensure specificity.
        attr_pattern (re.pattern): Compiled regex pattern to extract attribute
            names and values from the attribute column.
        fields (re.pattern): Compiled regex pattern to determine the number of
            fields (key, value pairs) in the attribute column.
        field_delimiter (re.pattern): Compiled regex pattern of the appropriate
            delimiter separating the key, value pairs.

    """
    col_pattern = None
    fields = None
    field_delimiter = None

    if file_type == 'gff3':
        # check if gff3 can really contain  "
        col_pattern = re.compile(r'([A-z0-9:;,+\-=._()" ]+)')
        fields = re.compile(r'([A-z0-9:,+\-=._()" ]+)')
        field_delimiter = re.compile(r'(=)')
    if file_type == 'gtf':
        col_pattern = re.compile(r'([A-z0-9:;+\-._()\" ]+)')
        fields = re.compile(r'([A-z0-9:,+\-._()\" ]+)')
        field_delimiter = re.compile(r'(;)')

    attr_pattern = re.compile(r'([^\=\;\"\s]+)')
    return col_pattern, attr_pattern, fields, field_delimiter


def valid_attributes(cols, field_pattern, field_delimiter, line_number):
    """Function that checks the validity of the gff/gtf attribute column.

    Args:
        cols (list): A gtf3 or gff line split into different columns.
        field_pattern (re.pattern): Compiled regex pattern returned from
            get_regex_pattern() to extract the fields (key, value pairs)
            in the attribute column.
        field_delimiter (re.pattern): Compiled regex pattern returned from
            get_regex_pattern() of the appropriate delimiter separating the
            key, value pairs.
        line_number (int): Number of the line to check for loggin purposes.

    Returns:
        check (boolean): Boolean indicating a failed or succcessfull check.

    """
    attribs = cols[ATTR_IDX]
    fields = field_pattern.findall(attribs)
    delims = field_delimiter.findall(attribs)
    check = len(fields) == len(delims)

    if not check:
        LOGGER.warning('Line: %d contains an invalid attribute column',
                       line_number)
    return check


def write_to_file(file_name, file_type, features):
    """Function that writes transcripts and exons to a gtf file.

     Args:
         file_name (str): The path of the file to be written.
         file_type (str): String idicating the file type. Has to be either
            'gff3' or 'gtf'.
         features (dict): A dict mapping the unique transcript id to objects of
             the feature class.

     """
    check_file_type_args(file_type)
    with open(file_name, 'w') as fout:
        for key, transcript in features.items():
            fout.write(transcript.get_as_string(file_type))


# TODO: Add docsting for the classes and methods
if sys.version_info[0] > 2:
    class Feature(ABC):
        """A class abstraction of a gff3/gtf genomic feature. Inherits from ABC.
        Implements most of the attributes of a genomic feature that are different
        for each type of feature. Class Transript and Exon inherit from it.

        Args:
            data (list): Columns of gff3/gtf line
            feature_id (str): A string holding a unique identifier.

        """
        def __init__(self, data, feature_id):
            # initilize variables and fill with some data
            self.name = feature_id
            self.score = data[SCORE_IDX]
            self.strand = data[STRAND_IDX]
            self.phase = data[FRAME_IDX]
            self.start = data[START_IDX]
            self.stop = data[STOP_IDX]
            self.attributes = data[ATTR_IDX]


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


        @abstractmethod
        def get_feature_id(self, attributes):
            """Getter method to return the unique feature id. Needs to be
            implemented by classes inherting from Feature.
            """
            pass
else:
    class Feature:
        """A class abstraction of a gff3/gtf genomic feature. Inherits from ABC.
        Implements most of the attributes of a genomic feature that are different
        for each type of feature. Class Transript and Exon inherit from it.

        Args:
            data (list): Columns of gff3/gtf line
            feature_id (str): A string holding a unique identifier.

        """
        __metaclass__ = ABCMeta
        def __init__(self, data, feature_id):
            # initilize variables and fill with some data
            self.name = feature_id
            self.score = data[SCORE_IDX]
            self.strand = data[STRAND_IDX]
            self.phase = data[FRAME_IDX]
            self.start = data[START_IDX]
            self.stop = data[STOP_IDX]
            self.attributes = data[ATTR_IDX]


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


        @abstractmethod
        def get_feature_id(self, attributes):
            """Getter method to return the unique feature id. Needs to be
            implemented by classes inherting from Feature.
            """
            pass



class Transcript(Feature):
    """A class abstraction of a gff3/gtf transcript . Inherits from Feature.
    Acts as container for data and Exon objects.

    Args:
        data (list): Columns of gff3/gtf line containing transcript
            information

    """
    def __init__(self, data):
        # initilize the Feature super class and inherit from it.
        super().__init__(data, self.get_feature_id(data[ATTR_IDX]))
        # initilize variables and fill with some data
        self.chrom = data[CHROM_IDX]
        self.source = data[SOURCE_IDX]
        self.exons = []
        self.exon_order = []
        self.number_exons = 0


    def __eq__(self, other):
        """Implentation of the comparison method. This method will compare the
        intron chain of two transcripts. It should be rather seen as a
        convienient way to compare the intron chain of two transcripts than a
        ground truth comparion on the whole transcript. It does not require
        exact matches of 5' and 3' ends of the fist and last exon. Returns true
        when both transcripts have no introns.

        Args:
            other (Transcript): Another transcript of the class Transcript.

        Retuns:
            check (boolean): Boolean indicating matching or not matching
                intron chains.

        """
        # transcripts can only be equal when located on the same strand and
        # chromosome and have an identical splice junctions/intron chain.
        if self.strand != other.strand:
            return False
        if self.chrom != other.chrom:
            return False
        if self.number_exons != other.number_exons:
            return False
        return self.compare_intron_chain(other)


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
            id_str_1 = ''.join(map(str, self.get_exon_start()[1:]))
            id_str_2 = ''.join(map(str, self.get_exon_end()[:-1]))
            return self.strand + self.chrom + id_str_1 + id_str_2
        # FIXME not an ideal solution. Need to come up with a way how to
        # deal with single exon transcripts
        return self.strand + self.chrom + '_1_exon'


    def compare_intron_chain(self, transcript, mode='full'):
        """Method to compare the intron chain of this transcript with that from
        another. This method only compares if the junctions share the same
        positions. It can either compare the full chain (mode='full') or
        check if another transcript is a super (mode='superset') or subset
        (mode='subset')of this object. It does not require exact matches of 5'
        and 3' ends of the fist and last exon. Returns true when both
        transcripts have no introns.

        Args:
            transcript (Transcript): An object of the class Transcript.
            mode (string): String indicating the type of check that should be
                performed. Needs to be either'full', 'superset' or 'subset'.
                Defaults to 'full'.

        Returns:
            check (boolean): Boolean indicating a failed or succcessfull check.

        Raises:
            ValueError: Raises an error when mode does not equal to either
                'full', 'superset' or 'subset'.

        """
        exons_self = set(self.get_exon_start()[1:])
        exons_self.update(self.get_exon_end()[:-1])

        exons_other = set(transcript.get_exon_start()[1:])
        exons_other.update(transcript.get_exon_end()[:-1])

        if mode == 'full':
            common = exons_self & exons_other
            n_junctions = max(len(exons_self), len(exons_other))
            return len(common) == n_junctions
        elif mode == 'subset':
            return exons_other.issubset(exons_self)
        elif mode == 'superset':
            return exons_other.issuperset(exons_self)
        else:
            raise ValueError('Invalid mode selected. Mode needs to be either'
                             'equal to \'full\' or \'subset\'.')


    def add_exon(self, data):
        """Setter method to add exons that belong to a transript.

        Args:
            data (list): Columns of gff3/gtf line containing exon
                information

        """
        # exon_number is needed in case no exon number is provided by the
        # annotation file
        exon_number = len(self.exons) + 1
        exon = Exon(data, exon_number)
        self.exons.append(exon)
        self.exon_order = self.sort_exons()
        self.number_exons = len(self.exons)

    def sort_exons(self):
        """Method to sort the exons according to ascending starting position.
        This method will be called after an exon has be added.

        """
        start_l = []
        for e in self.exons:
            start_l.append(e.start)
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

    def get_feature_id(self, attributes, verbose=False):
        """Getter method that returns a unique feature id contained in the
        attributes ('ID' or 'transcript_id').

        Returns:
            feture_id (str): The string of an unique id

        Raises:
            AttributeError: Will be raised when no 'ID' or 'transcript id'
                field is present in the attributes.

        """
        if verbose:
            print(attributes)
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

        for e in self.exons:
            line = ('%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s' %
                    (first_cols, 'exon', e.start, e.stop, e.score, e.strand,
                     e.phase, e.get_attributes(file_type)))
            lines.append(line)

        return '\n'.join(lines) + '\n'


class Exon(Feature):
    """A class abstraction of a gff3/gtf exon . Inherits from Feature.
    Acts as container for data.

    Args:
        data (list): Columns of gff3/gtf line containing transcript
            information
        exon_number (int): A number indicating the number of the exon.

    """
    def __init__(self, data, exon_number):
        self.parent = self.get_parent(data[ATTR_IDX])
        self.exon_number = exon_number
        super().__init__(data, self.get_feature_id(data[ATTR_IDX]))

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


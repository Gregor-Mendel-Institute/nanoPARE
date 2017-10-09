from collections import OrderedDict
# TODO needs testing and some reworking/bugfixing

def get_file_content(f):
    """Function that returns all lines contained in a the file 'f'.

    Args:
        f (str): A file path

    Returns:
        list: The lines of the file.

    """
    with open(f) as fin:
        return fin.readlines()


def is_feature(query, subject):
    """Funtion that checks if a query is contained in a string.

    Args:
        query (str): The substring that should be matched.
        subject (str): The GTF file line that should be matched against.

    Returns:
        bool: True if '\tquery\t' is contained in the subject string.

    """
    # add tabs to the query string to avoid feature column unspecific matches
    q = '\t' + query + '\t'
    return q in subject


def tokenize_string(s, delimiter):
    """Function that splits a string into substrings.

    Args:
        s (str): The string that should be split.
        delimiter (str): The delimiter that should be used for splitting.

    Returns:
        list: Returns a list of strings.

    """
    # remove trailing whitepsaces
    s = s.strip()
    # split and return
    return s.split(delimiter)


def process_attributes(attr):
    """Function that proccesses the attribute column of a GTF file and returns
    a dictionay with key: value mappings.

    Args:
        attr (str): The attribute column of a GTF file.

    Returns:
        dictionay: A dictionay mapping GTF attibutes to their values.

    """
    # replace " chars denoting the attribute values
    attr = attr.replace('"', '')
    # remove trailing ; to avoid an empty elements after string splitting
    if attr.endswith(";"):
        attr = attr[:-1]
    # split string
    tokens = tokenize_string(attr, ';')
    # remove trailing whitespaces
    tokens = list(map(str.strip, tokens))
    # split attributes at the whitespaces between key: values
    tokens = [x.split(' ') for x in tokens]
    # return as dictionary
    return OrderedDict(tokens)


def is_float(x):
    """Function that checks if a string is actually a floating point number.

    Args:
        attr (str): The string to test.

    Returns:
        boolean: True/False if float

    """
    try:
        a = float(x)
    except ValueError:
        return False
    else:
        return True

def is_int(x):
    """Function that checks if a string is actually an integer.

    Args:
        attr (str): The string to test.

    Returns:
        boolean: True/False if integer

    """
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a is b


def process_gft_line(line, id_field):
    # TODO needs more comments.
    """Function that processes one line of a GTF file, discards unnecessary
    information and returns a dictionary with key: value mappings.

    Args:
        line (str): The line to be processed (e.g. transcript,exon)
        feature_type (str): A string containing the feature type
        (e.g. exons, 5'UTR. etc)

    Returns:
        dictionary: A dictionary mapping the different attributes contained
                    in one line to their values.
    """
    header = ['chromosome',
              'source',
              'feature',
              'start',
              'end',
              'score',
              'strand',
              'frame',
              'attributes']

    tokens = tokenize_string(line, '\t')
    d = OrderedDict(zip(header, tokens))
    d.update(process_attributes(d['attributes']))
    del d['attributes']

    for x in d:
        t = d[x]
        if is_int(t):
            d[x] = int(t)
        if is_float(t):
            d[x] = float(t)

    return d[id_field], d['feature'], d


def parse_gtf(file_name, id_field, super_feature, sub_features, return_type,
              start_line=0, stop_line=0):
    """Function that parses a GTF file.

    Args:
        file_name (str): The path to the GTF file that should be parsed.
        id_field (str): The name of the unique identifier in the GTF
                         (e.g. transcript_id)
        super_feature (str): The name of the feature by which the subfeatures
                             should be grouped (e.g. transcript, mRNA, gene).
        sub_features (list[str]): A list of subfeatures that should be mapped
                                  to a super_feature (e.g. exon or 5'UTR).
        start_line (int): The line number from which on the file should be
                          parsed.
        stop_line (int): The line number until which the file should be parsed.

    Returns:
        dictionary: A dictionary mapping the id_field to a dictionary
                    containing attribute names and values.

    """
    # initialize some empty variables
    features = dict()
    feat = None
    # read in the file
    lines = get_file_content(file_name)
    # process the file line by line and add the features or subfeatures to a
    # dictionary
    for idx, line in enumerate(lines):
        # skipp lines till start_line
        if idx < start_line:
            continue
        next_line = idx + 1
        # check if the line is a superfeature and process it accordingly
        if is_feature(super_feature, line):
            # process the line
            unique_id, feat_type, data = process_gft_line(line, id_field)
            # initialize a Feature object and add the values
            feat = Feature(unique_id, feat_type, data)
        # check if the line is a subfeature and process it accordingly
        for sub_feature in sub_features:
            if is_feature(sub_feature, line):
                # process the line
                unique_id, feat_type, data = process_gft_line(line, id_field)
                # add subfeatures to the initialized feature object
                sub_feat = Feature(unique_id, feat_type, data)
                feat.add_subfeature(sub_feat)

        # add the feature to dictionary mapping id_field: {features} when the
        # next line is the beginning of a new feature, the file has ended or
        # the maximum number of lines to process has been reached
        if (next_line in {len(lines), stop_line} or
                is_feature(super_feature, lines[next_line]) and
                feat is not None):
            #print(feat.name)
            features[feat.name] = feat
            # reset variable to allow initializing of a new Feature object in
            # the next iteration
            feat = None
    return features


def get_file_content(f):
    """Function that returns all lines contained in a the file 'f'.

    Args:
        f (str): A file path

    Returns:
        list: The lines of the file.

    """
    with open(f) as fin:
        return fin.readlines()


def write_gtf(file_name, features):

    with open(file_name, 'w') as fout:
        for key in features:
            feat = features[key]
            write_line(fout, feat.data)
            for sfeat in feat.subfeatures:
                write_line(fout, sfeat.data)

def write_line(f, d):
    header = ['chromosome',
              'source',
              'feature',
              'start',
              'end',
              'score',
              'strand',
              'frame']

    attributes = ''
    cols = []
    for key in d:
        value = d[key]
        if key in header:
            cols.append(str(value))
        else:
            attributes = attributes + key + ' "' + str(value) + '"' + '; '
    line = '\t'.join(cols) + '\t' + attributes[:-2] + '\n'
    f.write(line)


class Feature(object):
    "Generic Feature node."
    def __init__(self, name, feat_type, data, subfeatures=None):
        self.name = name
        self.feat_type = feat_type
        self.data = data
        self.subfeatures = []
        if subfeatures is not None:
            for s in subfeatures:
                self.add_subfeature(s)

    def __repr__(self):
        return (self.name + '_' + self.feat_type)

    def add_data(data):
        self.data = data

    def add_subfeature(self, feat):
        assert isinstance(feat, Feature)
        self.subfeatures.append(feat)
'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Blah XXX 
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import csv
from intervaltree import Interval, IntervalTree


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_CNV_FILE_ERROR = 3
DEFAULT_VERBOSE = False
DEFAULT_BND_WINDOW = 50
DEFAULT_CONF_STEPS = 10
PROGRAM_NAME = "repagree"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Blah XXX'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('--overlap',
                        metavar='OVERLAP',
                        type=float,
                        help='proportion of overlap required between two CNVs to be considered equal')
    parser.add_argument('--steps',
                        metavar='STEPS',
                        type=int,
                        default=DEFAULT_CONF_STEPS,
                        help='number of steps in confidence threshold')
    parser.add_argument('cnv_file',
                        metavar='CNF_FILE',
                        type=str,
                        help='Input CNV file')
    return parser.parse_args()


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))

def read_variants(cnv_filename):
    sample_ids = set()
    conf_scores = []
    variants = []
    with open(cnv_filename) as file:
        reader = csv.DictReader(file)
        for row in reader:
            conf_scores.append(float(row['conf']))
            variants.append(row)
            sample_ids.add(row['sampleID'])
    if len(sample_ids) == 2:
        sample_list = list(sample_ids)
        min_conf = min(conf_scores)
        max_conf = max(conf_scores)
        return sample_list[0], sample_list[1], min_conf, max_conf, variants
    else:
        exit_with_error("Wrong number of sample IDs: {}".format(sample_ids), EXIT_CNV_FILE_ERROR)

class Intervals(object):
    def __init__(self):
        self.chroms = {}

    def insert(self, chrom, start, end, val):
        if chrom not in self.chroms:
            self.chroms[chrom] = IntervalTree()
        self.chroms[chrom][start:end] = val 

    def lookup(self, chrom, start, end):
        if chrom in self.chroms:
            return self.chroms[chrom][start:end]
        else:
            return set()

def populate_intervals(sample_id, variants):
    intervals = Intervals()
    for row in variants:
        if row['sampleID'] == sample_id:
            val = (row['cn'], float(row['conf']))
            intervals.insert(row['chr'], int(row['start']), int(row['end']), val)
    return intervals


def overlap_filter(start1, end1, start2, end2, min_overlap):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start < overlap_end:
        overlap_size = float((overlap_end - overlap_start) + 1)
        cnv1_size = (end1 - start1) + 1
        cnv2_size = (end2 - start2) + 1
        cnv1_overlap = overlap_size / cnv1_size
        cnv2_overlap = overlap_size / cnv2_size
        return cnv1_overlap >= min_overlap and cnv2_overlap >= min_overlap
    return False


def filter_matches(overlap_threshold, cn, start, end, conf_threshold, variants):
    result = 0
    for variant in variants:
        this_cn, this_conf = variant.data
        if this_cn == cn and this_conf >= conf_threshold and overlap_filter(start, end, variant.begin, variant.end, overlap_threshold):
            result += 1
    return result

def agreement(overlap_threshold, step_size, max_conf, sample_id, intervals, variants):
    result = []
    for conf_threshold in range(0, max_conf + step_size, step_size):
        num_agree_cnvs = 0
        total_variants_above_threshold = 0
        for row in variants:
            this_conf = float(row['conf'])
            if this_conf >= conf_threshold:
                total_variants_above_threshold += 1
                if row['sampleID'] == sample_id:
                    this_cn = row['cn']
                    this_start = int(row['start'])
                    this_end = int(row['end'])
                    overlaps = intervals.lookup(row['chr'], this_start, this_end)
                    num_agree_cnvs += filter_matches(overlap_threshold, this_cn, this_start, this_end, conf_threshold, overlaps)
        denominator = total_variants_above_threshold - num_agree_cnvs
        if denominator > 0:
            proportion = num_agree_cnvs / denominator 
            result.append((conf_threshold, proportion))
    return result


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    sample_a, sample_b, min_conf, max_conf, variants = read_variants(options.cnv_file)
    step_size = int(max_conf / options.steps)
    interval_tree = populate_intervals(sample_b, variants)
    points = agreement(options.overlap, step_size, int(max_conf + 1), sample_a, interval_tree, variants)
    print(points)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()

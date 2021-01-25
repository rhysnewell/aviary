#!/usr/bin/env python
###############################################################################
# 
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
__author__ = "Rhys Newell"
__copyright__ = "Copyright 2020"
__credits__ = ["Rhys Newell"]
__version__ = "1.0.0"
__license__ = "GPL3"
__maintainer__ = "Rhys Newell"
__email__ = "rhys.newell near hdr.qut.edu.au"
__status__ = "Development"

###############################################################################
# System imports
import sys
import argparse
import logging
import os
import datetime
from pathlib import Path

# Function imports
import numpy as np
import pandas as pd
from Bio import SeqIO

# Debug
debug = {
    1: logging.CRITICAL,
    2: logging.ERROR,
    3: logging.WARNING,
    4: logging.INFO,
    5: logging.DEBUG
}

###############################################################################
############################### - Exceptions - ################################


class BadTreeFileException(Exception):
    pass


###############################################################################                                                                                                                      [44/1010]
################################ - Functions - ################################
def main():
    ############################ ~ Main Parser ~ ##############################
    vamb_options = argparse.ArgumentParser(
        'write_vamb_bins',
        description='Bin out the results of vamb',
        formatter_class=CustomHelpFormatter,
        epilog='''
                                        ~ write_vamb_bins ~
            How to use write_vamb_bins:

            write_vamb_bins --reference assembly.fasta --clusters vamb_clusters.tsv

            NOTE: Provided defaults are for use with aviary and should be changed for your specific use

            ''')
    vamb_options.add_argument('--version',
                             action='version',
                             version=__version__,
                             help='Show version information.')
    vamb_options.add_argument(
        '--verbosity',
        help=
        '1 = critical, 2 = error, 3 = warning, 4 = info, 5 = debug. Default = 4 (logging)',
        type=int,
        default=4)
    vamb_options.add_argument('--log',
                             help='Output logging information to file',
                             default=False)

    ########################## ~ sub-parser ~ #########################

    vamb_options.add_argument('--reference',
                              help='The assembly file to be binned',
                              dest='assembly',
                              default="data/vamb_bams/renamed_contigs.fasta")

    vamb_options.add_argument('--clusters',
                              help='The vamb clusters',
                              dest='clusters',
                              default="data/vamb_bins/clusters.tsv")

    vamb_options.add_argument('--min_size',
                              help='Minimum bin size',
                              dest='min_size',
                              default=200000)

    vamb_options.add_argument('--output',
                              help='The output directory',
                              dest='output',
                              default='data/vamb_bins/')

    vamb_options.set_defaults(func=vamb)
    ###########################################################################
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parsing input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    args = vamb_options.parse_args()
    time = datetime.datetime.now().strftime('%H:%M:%S %d-%m-%Y')

    if args.log:
        if os.path.isfile(args.log):
            raise Exception("File %s exists" % args.log)
        logging.basicConfig(
            filename=args.log,
            level=debug[args.verbosity],
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(
            level=debug[args.verbosity],
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info("Time - %s" % (time))
    logging.info("Command - %s" % ' '.join(sys.argv))

    args.func(args)


def vamb(args):
    min_bin_size = int(args.min_size)
    prefix = args.output
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    bins = {}
    with open(args.clusters, 'r') as vamb_file:
        for line in vamb_file:
            line = line.split()
            try:
                bins[line[0]].append(line[1])
            except KeyError:
                bins[line[0]] = [line[1]]

    assembly = SeqIO.to_dict(SeqIO.parse(args.assembly, "fasta"))

    logging.info("Writing bins...")
    max_cluster_id = max(bins.keys())
    for (bin, contigs) in bins.items():
        if bin != -1:
            # Calculate total bin size and check if it is larger than min_bin_size
            bin_length = sum([len(assembly[contig].seq) for contig in contigs])
            if bin_length >= min_bin_size:
                with open(prefix + '/vamb_bin.' + str(bin) + '.fna', 'w') as f:
                    for contig in contigs:
                        write_contig(contig, assembly, f)

        else:
            # Get final bin value
            max_cluster_id += 1
            # Rescue any large unbinned contigs and put them in their own cluster
            for contig in contigs:
                if len(assembly[contig].seq) >= min_bin_size:
                    with open(prefix + '/vamb_bin.' + str(max_cluster_id) + '.fna', 'w') as f:
                        write_contig(contig, assembly, f)

    Path('data/vamb_bins/done').touch()


def write_contig(contig, assembly, f):
    seq = assembly[contig]
    fasta = ">" + seq.id + '\n'
    fasta += str(seq.seq) + '\n'
    f.write(fasta)

###############################################################################
################################ - Classes - ##################################


class CustomHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if '%(default)' not in action.help:
            if action.default != '' and \
               action.default != [] and \
               action.default != None \
               and action.default != False:
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [
                        argparse.OPTIONAL, argparse.ZERO_OR_MORE
                    ]

                    if action.option_strings or action.nargs in defaulting_nargs:

                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
        return h

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])



if __name__ == '__main__':

   sys.exit(main())

#!/usr/bin/env python
# Import the necessary modules
# OS is used for file/folder manipulations
import os
# Subprocess->call is used for making system calls
import subprocess
# Errno is used in the file creation command  - I think it's similar to the $! variable in Perl
import errno
# Glob finds all the path names matching a specified pattern according to the rules used by the Unix shell
from glob import glob
# Shutil is useful for file moving/copying
import shutil
# System tools
import sys
# Can import date, and calculate length of run, etc.
import time
# Multiprocessing module
from multiprocessing import Pool
# Regex module
import re
# JSON module for reading and printing variables in .json format
# import json
# Sequence parsing module from BioPython
from Bio import SeqIO
# Default dictionaries allow nested dictionaries
from collections import defaultdict
# Import custom modules
import SMALTcombined
import SMALT
import bamProcessorCombined
import bamProcessor
import rawMLST
import bamPysamStats
import bamPysamStatsCombined
import fastqCreator

__author__ = 'adamkoziol'


class GeneSippr(object):

    def sippr(self):
        pass

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """

        :param args:
        :param pipelinecommit:
        :param startingtime:
        :param scriptpath:
        """
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args.path, '')
        self.commit = str(pipelinecommit)
        self.start = startingtime
        self.homepath = scriptpath


if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Parser for arguments
    parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
    parser.add_argument('-v', '--version',
                        version='%(prog)s commit {}'.format(commit))
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        help='Path of .fastq(.gz) files to process. If not provided, the default path of '
                             '"path/sequences" will be used')
    parser.add_argument('-t', '--targetpath',
                        help='Path of target files to process. If not provided, the default path of '
                             '"path/targets" will be used')
    parser.add_argument('-m', '--miSeqPath',
                        help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', '--miseqfolder',
                        help='Name of the folder containing MiSeq run data')
    parser.add_argument('-r1', '--readLengthForward',
                        help='Length of forward reads to use. Can specify "full" to take the full length of '
                             'forward reads specified on the SampleSheet')
    parser.add_argument('-r2', '--readLengthReverse',
                        default=0,
                        help='Length of reverse reads to use. Can specify "full" to take the full length of '
                             'reverse reads specified on the SampleSheet')
    parser.add_argument('-c', '--customSampleSheet',
                        help='Path of folder containing a custom sample sheet (still must be named "SampleSheet.csv")')
    parser.add_argument('-P', '--projectName',
                        help='A name for the analyses. If nothing is provided, then the "Sample_Project" field '
                             'in the provided sample sheet will be used. Please note that bcl2fastq creates '
                             'subfolders using the project name, so if multiple names are provided, the results '
                             'will be split as into multiple projects')
    parser.add_argument('-16S', '--16Styping',
                        action='store_true',
                        help='Perform 16S typing. Note that for analyses such as MLST, pathotyping, '
                             'serotyping, and virulence typing that require the genus of a strain to proceed, '
                             '16S typing will still be performed')
    parser.add_argument('-M', '--Mlst',
                        action='store_true',
                        help='Perform MLST analyses')
    parser.add_argument('-Y', '--pathotYping',
                        action='store_true',
                        help='Perform pathotyping analyses')
    parser.add_argument('-S', '--Serotyping',
                        action='store_true',
                        help='Perform serotyping analyses')
    parser.add_argument('-V',
                        '--Virulencetyping',
                        action='store_true',
                        help='Perform virulence typing analyses')
    parser.add_argument('-a', '--armi',
                        action='store_true',
                        help='Perform ARMI antimicrobial typing analyses')
    parser.add_argument('-r', '--rmlst',
                        action='store_true',
                        help='Perform rMLST analyses')
    parser.add_argument('-d', '--detailedReports',
                        action='store_true',
                        help='Provide detailed reports with percent identity and depth of coverage values '
                             'rather than just "+" for positive results')
    parser.add_argument('-C', '--customTargetPath',
                        help='Provide the path for a folder of custom targets .fasta format')

    # TODO Add custom cutoffs
    # TODO Assert .fastq files present in provided folder
    # TODO Don't touch .fastq(.gz) files

    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    GeneSippr(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'
    # print json.dumps(seqdict, sort_keys=True, indent=4, separators=(',', ': '))

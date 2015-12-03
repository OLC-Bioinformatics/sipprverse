__author__ = 'akoziol'

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
# Regex
import re
# System tools
import sys
# Can import date, and calculate length of run, etc.
import time
# Multiprocessing module
from multiprocessing import Pool
# Numerical python - used in the parsing of vcf files
import numpy
# Math module - used in the parsing of vcf files
import math
# JSON module for reading and printing variables in .json format
import json
# Default dictionaries allow nested dictionaries
from collections import defaultdict
# Argument parser for user-inputted values, and a nifty help menu
from argparse import ArgumentParser

#Parser for arguments
parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
parser.add_argument('-v', '--version', action='version', version='%(prog)s v1.0')
parser.add_argument('-p', '--path', required=True, help='Specify input directory')
parser.add_argument('-s', '--sequencePath', required=False, help='Path of .fastq(.gz) files to process. If not '
                    'provided, the default path of "path/sequences" will be used')
parser.add_argument('-t', '--targetPath', required=False, help='Path of target files to process. If not '
                    'provided, the default path of "path/targets" will be used')
# parser.add_argument('-l', '--readLength', required=True, help='Specify list of read lengths to be used e.g. 18, 19, 20, 21, 22')
# parser.add_argument('-f', '--foldCoverage', required=True, help='Specify list of fold coverage values to be used e.g. 1, 1.5, 2, 2.5, 5')
# parser.add_argument('-k', '--kmerLength', required=True, help='Specify list of kmer lengths to be used e.g. 5, 7, 11')

# Get the arguments into a list
args = vars(parser.parse_args())

# Define variables from the arguments - there may be a more streamlined way to do this
path = os.path.join(args['path'], "")
# If these variables are not defined, set to default values
if args['sequencePath']:
    sequencePath = os.path.join(args['sequencePath'], "")
else:
    sequencePath = "%ssequences/" % path

if args['targetPath']:
    targetPath = os.path.join(args['targetPath'], "")
else:
    targetPath = "%stargets/" % path

targets = glob("%s*.fasta" % targetPath)

for target in targets:
    print os.path.basename(target).split(".")[0]


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def indexTargetsProcesses(targets):
    """Faidx multiprocessing helper function"""
    print '\nProcessing targets with faidx'
    # Initialise the args list
    indexArgs = []
    if __name__ == '__main__':
        # Initialise the pool of processes - it defaults to the number of processors
        indexPool = Pool()
        # for target in targets:
            # indexArgs.append(target)
        indexPool.map(indexTargets, targets)


def indexTargets((target)):
    global targetPath
    print target
    # Format the target names properly
    filename = os.path.split(target)[1]
    fileNoExt = filename.split(".")[0]
    # Create a new path to be created (if necessary) for the generation of the range of k-mers
    # indexPath = "%stargets" % path
    # Call the make_path function to make folders as necessary
    # make_path(indexPath)
    indexFileSMI = "%s.smi" % fileNoExt
    indexFileSMA = "%s.sma" % fileNoExt
    print filename, fileNoExt, indexFileSMI
    # Index the appropriate files
    # if not os.path.isfile("%s/%s" % (targetPath, indexFileSMI)):
    indexCommand = "smalt index -k 20 -s 1 %s %s" % (fileNoExt, target)
    subprocess.call(indexCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    try:
        shutil.move(indexFileSMI, targetPath)
    except OSError, shutil.Error:
        pass
    try:
        shutil.move(indexFileSMA, targetPath)
    except OSError, shutil.Error:
        pass
    sys.stdout.write('.')

    # else:
    #     sys.stdout.write('.')

indexTargetsProcesses(targets)
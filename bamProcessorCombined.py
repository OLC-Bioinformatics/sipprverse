import subprocess
import os
import sys
import time
import errno
# Multiprocessing module
from multiprocessing import Pool
from glob import glob
from collections import defaultdict

__author__ = 'akoziol'

# Initialise the count used in the dotter function
count = 0


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    # Use a global variable
    global count
    if count <= 80:
        sys.stdout.write('.')
        count += 1
    else:
        sys.stdout.write('\n[%s] .' % (time.strftime("%H:%M:%S")))
        count = 1


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL.
    :param inpath: string of the supplied path
    """
    try:
        os.makedirs(inpath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)


def sortingprocesses(seqdict, analysistype):
    """
    Multiprocessing for sorting bam files
    :param seqdict: dictionary containing important file and path information
    :param analysistype: string of the current analysis type
    """
    # Initialise the arguments and Pool variables
    sortingprocessesargs = []
    sortingprocessespool = Pool()
    # Iterate through the strains
    for strain in seqdict:
        # Pull the bait type from seqdict
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        # Retrieve name and path of the concatenated file to use in smalt mapping
        catfile = seqdict[strain]["concatenatedTargets"][analysistype]
        # Find the location of the fastq files
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        # Append a tuple of the arguments to the arguments list
        sortingprocessesargs.append((catfile, fastqdir))
    # Map the arguments to the sorting function
    sortingprocessespool.map(sorting, sortingprocessesargs)


def sorting((target, fastqdir)):
    """Performs samtools sort to return a sorted bam file"""
    # Strip off the path and the file extension from the targets
    targetname = os.path.basename(target).split(".")[0]
    # Define the folder to find inputs and store outputs
    outdir = "%s/%s" % (fastqdir, targetname)
    # Find bam files in the input folder
    bamfiles = glob("%s/*.bam" % outdir)
    # Iterate through the files
    for bamfile in bamfiles:
        # Ignore any files that have been sorted previously
        if "_sorted" not in bamfile:
            # Define the name for the sorted files
            sortedname = "%s_sorted" % bamfile.split(".")[0]
            # If sorting has not already occurred
            if not os.path.isfile("%s.bam" % sortedname):
                # Define the samtools sort call
                bamsort = "samtools sort %s %s" % (bamfile, sortedname)
                # Set /dev/null
                fnull = open(os.devnull, 'wb')
                # Run the call
                subprocess.call(bamsort, shell=True, stdout=fnull, stderr=fnull)
        dotter()


def bamindexingprocesses(seqdict, analysistype):
    """
    Multiprocessing for indexing bam files
    :param seqdict: dictionary containing important file and path information
    :param analysistype: string of the current analysis type
    """
    # Initialise argument and Pool variables
    bamindexingargs = []
    bamindexingpool = Pool()
    for strain in seqdict:
        # Determine bait type for this analysis type
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        # Retrieve name and path of the concatenated file to use in smalt mapping
        catfile = seqdict[strain]["concatenatedTargets"][analysistype]
        # Find the location of the fastq files
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        # Append a tuple of the arguments to the arguments list
        bamindexingargs.append((catfile, fastqdir))
    # Map the arguments list to the indexing function
    bamindexingpool.map(bamindexing, bamindexingargs)


def bamindexing((target, fastqdir)):
    """Indexes the sorted bam files in order to visualize the assemblies with tablet - note this is OPTIONAL"""
    """Indexes the sorted bam files"""
    # Extract the name of each target from the supplied target variable
    targetname = os.path.basename(target).split(".")[0]
    # Set the input/output directory
    outdir = "%s/%s" % (fastqdir, targetname)
    # Find the sorted bam files
    bamfiles = glob("%s/*_sorted.bam" % outdir)
    for bamfile in bamfiles:
        # Set the name of output file
        baifile = bamfile.split(".")[0] + ".bai"
        # Only run the indexing if the .bai file is not present in the output folder
        if not os.path.isfile(baifile):
            # Define the indexing call
            bamindex = "samtools index %s %s" % (bamfile, baifile)
            # Define /dev/null
            fnull = open(os.devnull, 'wb')
            # Run the call
            subprocess.call(bamindex, shell=True, stdout=fnull, stderr=fnull)
        dotter()

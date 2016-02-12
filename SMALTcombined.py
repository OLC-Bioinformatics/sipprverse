import subprocess
import os
import sys
import time
import errno
# Multiprocessing module
from multiprocessing import Pool

__author__ = 'adamkoziol'

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


def smaltindextargets(target, targetpath):
    """
    Function that indexes target files with SMALT
    :param target: string of concatenated target file
    :param targetpath: string of path of the target files
    """
    # Format the target names properly
    filename = os.path.split(target)[1]
    filenoext = filename.split(".")[0]
    # Index the appropriate files
    if not os.path.isfile("%s/%s.smi" % (targetpath, filenoext)):
        # Define the indexing command
        indexcommand = "cd %s && smalt index -k 20 -s 1 %s %s" % (targetpath, filenoext, target)
        # print indexcommand
        # Define /dev/null
        fnull = open(os.devnull, 'wb')
        # Run the command
        subprocess.call(indexcommand, shell=True, stdout=fnull, stderr=fnull)
    dotter()


def smaltmappingprocesses(seqdict, analysistype, mappername):
    """SMALT reference mapper helper function
    :param seqdict: dictionary containing important file and path information
    :param analysistype: string of the current analysis type
    :param mappername: string of the current reference mapping software
    """
    mappingprocessesargs = []
    mappingprocessespool = Pool()
    # Iterate through all the strains
    for strain in seqdict:
        # Get the value for baittype for each strain stored in seqdict
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        # Retrieve name and path of the concatenated file to use in smalt mapping
        catfile = seqdict[strain]["concatenatedTargets"][analysistype]
        # Find the list of baited fastq files to process
        baitedfastqlist = seqdict[strain]["bait"]["fastqFiles"][baittype]
        mappingprocessesargs.append((strain, catfile, mappername, baitedfastqlist))
    # Map the arguments to the reference mapping function
    mappingprocessespool.map(mapping, mappingprocessesargs)


def mapping((strain, target, mappername, baitedfastqlist)):
    """Performs the mapping of the reads to the targets"""
    # Initialise the smaltmap variable
    smaltmap = ""
    # Manipulate the target variable to give desired results
    targetnoext = target.split(".")[0]
    targetname = os.path.basename(target).split(".")[0]
    # Define the prefix to be appended to output bam files
    bamname = "%s_%s_%s" % (strain, targetname, mappername)
    # As the script allows paired end, or unpaired end fastq reads, both possibilities must be accounted for
    if len(baitedfastqlist) == 2:
        # Set the names of the forward and reverse fastq files from the supplied list
        forwardfastq, reversefastq = sorted(baitedfastqlist)
        # Determine the path of the fastq files by splitting off the file name
        fastqpath = os.path.split(forwardfastq)[0]
        # Set the path to store the outputs of the reference mapping
        outpath = "%s/%s" % (fastqpath, targetname)
        # Create the path if necessary
        make_path(outpath)
        # Perform the reference mapping if the analysis has not previously been performed
        # Note that this only checks to see if the desired file exists, it does not check if it size 0
        if not os.path.isfile("%s/%s.bam" % (outpath, bamname)):
            # Define the mapping call
            # The -d 0 flag allows all alignments or alignment pairings with the best score to be printed
            # The -S subst=-5 flag increases substitution penalties to -5 (up from a default of -2)
            smaltmap = "smalt map -o %s/%s.bam -f bam -n 24 -d 0 -l pe -S subst=-5 -x %s %s %s" \
                       % (outpath, bamname, targetnoext, forwardfastq, reversefastq)
    else:
        # Same as above, but with only one fastq file
        forwardfastq = baitedfastqlist[0]
        fastqpath = os.path.split(forwardfastq)[0]
        outpath = "%s/%s" % (fastqpath, targetname)
        make_path(outpath)
        if not os.path.isfile("%s/%s.bam" % (outpath, bamname)):
            # Define the mapping call
            # -d -1 is used with unpaired reads: if -d is set to a value < 0, all alignments are reported
            smaltmap = "smalt map -o %s/%s.bam -f bam -n 24 -d -1 -S subst=-5 -x %s %s" \
                       % (outpath, bamname, targetnoext, forwardfastq)
    # Define /dev/null
    fnull = open(os.devnull, 'wb')
    # Run the system call
    subprocess.call(smaltmap, shell=True, stdout=fnull, stderr=fnull)
    dotter()

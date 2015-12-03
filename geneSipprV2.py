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


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the metadata dictionaries
metadata = defaultdict(make_dict)
targetDict = {}


def loadTargetData():
    """Loads the dictionary of target data (if it exists)"""
    global path, targetDict, targetPath
    make_path("%sjsonFiles" % path)
    targetJSONfile = "%sjsonFiles/targetDict.json" % path
    if os.path.isfile(targetJSONfile):
        with open("%sjsonFiles/targetDict.json" % path) as targetJSON:
            targetDict = json.load(targetJSON)


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def timePrint(string):
    """Prints a formatted string with [hour:minute:second] text"""
    print "[%s] %s" % (time.strftime("%H:%M:%S"), string)

# Initialise count to 0
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


def foldererPrepProcesses(sampleList):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    foldererPrepArgs = []
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == '__main__':
        createfoldererPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleList:
            foldererPrepArgs.append(name)
        # This map function allows for multi-processing
        createfoldererPool.map(folderer, foldererPrepArgs)


def folderer(gzFile):
    """Uses gzip to decompress .gz.fastq files"""
    # Gzip each file
    gzipCommand = "gzip -d --force %s" % gzFile
    # System call
    subprocess.call(gzipCommand, shell=True,
                                 stdout=open(os.devnull, 'wb'),
                                 stderr=open(os.devnull, 'wb'))
    dotter()


def fileR(fileList):
    """Helper script that creates a set of the stain names created by stripping off parts of the filename.
    Hopefully handles different naming conventions (e.g. 2015-SEQ-001_S1_L001_R1_001.fastq(.gz),
    2015-SEQ-001_R1_001.fastq.gz, 2015-SEQ-001_R1.fastq.gz, 2015-SEQ-001_1.fastq.gz, and 2015-SEQ-001_1.fastq.gz
    all become 2015-SEQ-001)"""
    # Initialise the set
    fileSet = set()
    for seqPath in fileList:
        seqFile = os.path.split(seqPath)[1]
        # Search for the conventional motifs present following strain names
        # _S\d+_L001_R\d_001.fastq(.gz) is a typical unprocessed Illumina fastq file
        if re.search("_S\d+_L001", seqFile):
            fileSet.add(re.split("_S\d+_L001", seqFile)[0])
        # Files with _R\d_001.fastq(.gz) are created in the SPAdes assembly pipeline
        elif re.search("_R\d_001", seqFile):
            fileSet.add(re.split("_R\d_001", seqFile)[0])
        # _R\d.fastq(.gz) represents a simple naming scheme for paired end reads
        elif re.search("R\d.fastq", seqFile):
            fileSet.add(re.split("_R\d.fastq", seqFile)[0])
        # _\d.fastq(.gz) is always possible
        elif re.search("[-_]\d.fastq", seqFile):
            fileSet.add(re.split("[-_]\d.fastq", seqFile)[0])
        # .fastq is the last option
        else:
            fileSet.add(re.split(".fastq", seqFile)[0])
        dotter()
    # Convert to a list
    fileList = list(fileSet)
    return fileList


def fileGrabber(targetPath, sequencePath):
    """Prepares lists of files for use in subsequent functions"""
    # Target files are in path/targets
    # Must have a .*fa* extension
    timePrint("Acquiring target files")
    targetFiles = glob("%s*fa*" % targetPath)
    # Sequence files are in path/sequences
    # Check for .gz compressed files
    gzFiles = glob("%s*.gz" % sequencePath)
    if gzFiles:
        # Decompress the files
        timePrint("Decompressing files")
        foldererPrepProcesses(gzFiles)
        sys.stdout.write("\n")
    # Get all the .fastq files
    timePrint("Acquiring sequence files")
    sequenceFiles = glob("%s*.fastq" % sequencePath)
    # Get the unique file names
    timePrint("Acquiring sequence names")
    # Get unique names of sequences stripped from the rest of the file name
    sequenceNames = fileR(sequenceFiles)
    return targetFiles, sequenceFiles, sequenceNames


def indexTargets(targets):
    """Performs smalt index on the targets using the range of k-mers stored in the variable kmer"""
    global targetPath, path, targetDict
    timePrint('\nIndexing targets')
    # As far as I can tell, smalt index
    for target in targets:
        # Because there is a trailing slash on targetPath, shutil.move won't move files properly
        indexpath = "%stargets" % path
        # Format the target names properly
        filename = os.path.split(target)[1]
        fileNoExt = filename.split(".")[0]
        # Define the files created by the system call
        indexFileSMI = "%s/%s.smi" % (indexpath, fileNoExt)
        # Index the appropriate files only if necessary
        if not targetDict[fileNoExt] or not os.path.isfile(indexFileSMI):
            indexCommand = "smalt index -k 20 -s 1 %s %s" % (fileNoExt, target)
            # Perform the system call
            subprocess.call(indexCommand, shell=True,
                                          stdout=open(os.devnull, 'wb'),
                                          stderr=open(os.devnull, 'wb'))
            #  As far as I can tell, there is no option for smalt index to write files to any directory except the cwd
            #  Need to move the files to the appropriate path
            shutil.move("%s.smi" % fileNoExt, indexpath)
            shutil.move("%s.sma" % fileNoExt, indexpath)
            # Capture the indexing command for each target
            targetDict[fileNoExt] = indexCommand
            # Write the dictionary to file
            with open("%sjsonFiles/targetDict.json" % path, "wb") as targetJSONfile:
                # Print the JSON data to file
                output = json.dumps(targetDict, sort_keys=True, indent=4, separators=(',', ': '))
                targetJSONfile.write(output)
            dotter()
        else:
            dotter()


# def createOutputFiles():
#     """Parses the vcf files created above to create a handy summary table of mapping stats"""
#     print "\nCreating outputs"
#     make_path(outPath)
#     os.chdir(outPath)
#     outFile = open("SipprModelling_%s.csv" % start, "wb")
#     outFile.write("readLength\tfoldCoverage\ttarget\tkmerLength\tMedianQualityScore\t"
#                   "QualityScoreSD\tMedianFoldCoverage\tFoldCoverageSD\tMedianPercentID\tqualityMetric\n")
#     for rLength in readLength:
#             for fCov in foldCoverage:
#                 for target in targets:
#                     for size in kmer:
#                         total1 = 0
#                         sys.stdout.write('.')
#                         filename = os.path.split(target)[1]
#                         fileNoExt = filename.split(".")[0]
#                         megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, fileNoExt, size)
#                         filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
#                         vcfFile = megaName + "_sorted.vcf"
#                         newPath = "%s/%s" % (filePath, megaName)
#                         outputFile = "%s/%s" % (newPath, vcfFile)
#                         # Initialise the counter, which will be used to track lines in the vcf file - if positions in the
#                         # target are not mapped, then the position field will jump ahead of the counter
#                         count = 1
#                         # Initialise the arrays, which will keep track of the appropriate values for each dataset
#                         arrQual = []
#                         arrCov = []
#                         arrSum = []
#                         output = open(outputFile, "r")
#                         for line in output:
#                             # vcf files have 36 commented out lines at the top of each file - these are not necessary
#                             if re.search('#', line):
#                                 pass
#                             else:
#                                 total1 += 1
#                                 # Format of file
#                                 # CHROM	    POS	ID	REF	ALT	QUAL FILTER	INFO	                                   FORMAT
#                                 # adk-12	8	.	G	.	32.7	.	DP=1;AF1=0;AC1=0;DP4=0,1,0,0;MQ=29;FQ=-30.3	PL	0
#                                 # data[0] [1]  [2] [3]  [4] [5]    [6]  [7]
#                                 data = line.split("\t")
#                                 #target = data[0]
#                                 pos = data[1]
#                                 refSeq = data[3]
#                                 mapSeq = data[4]
#                                 qual = data[5]
#                                 # Depth of coverage is reported prior to the first ";"
#                                 dpLine = data[7].split(";")[0]
#                                 # For now, I'm skipping lines that indicated the presence of a possible indel
#                                 # - I may return to this later
#                                 if re.search("INDEL", dpLine):
#                                     pass
#                                 else:
#                                     # If the called base (mapSeq) is identical to the reference base (refSeq)
#                                     # - denoted by a ".", then set seq to equal refSeq, otherwise, pull the
#                                     # value of mapSeq for seq
#                                     avgQual = sum(arrQual)/total1
#                                     if mapSeq == ".":
#                                         seq = refSeq
#                                         match = 1
#                                     # This section corrects for the fact that during the conversion of bam files to vcf
#                                     # files, SNP calls and ambiguous calls look identical, except for the fact that for
#                                     # SNPs, the qualityScore (qual) tends to be higher than the surrounding bases,
#                                     # while ambiguous calls have a lower qualityScore - this loop screens for quality
#                                     # scores that are at least 10 lower than the score of the previous base
#                                     else:
#                                         if float(arrQual[-1] - 10) >= 0:
#                                             prevValue = float(arrQual[-1] - 10)
#                                         else:
#                                             prevValue = 0
#                                         if float(qual) <= prevValue:
#                                             seq = refSeq
#                                             match = 1
#                                         else:
#                                             # This attempts to catch if there are two ambiguous bases in a row;
#                                             # they will hopefully have the same value
#                                             if float(qual) == prevValue:
#                                                 seq = refSeq
#                                                 match = 1
#                                             else:
#                                                 # "True" SNPs seem to have increased qualityScore compared to the
#                                                 # surrounding values, this will catch that
#                                                 if float(qual) > prevValue:
#                                                     seq = mapSeq
#                                                     match = 0
#                                     # Strip the "DP=" from dpLine
#                                     DP = dpLine.split("=")[1]
#                                     #vcfData[pos] = (fileName, target, refSeq, mapSeq, DP)
#                                     # If pos > count, then there is a gap in the mapping (or a deletion, but ignoring
#                                     # this possibility for now). For my purposes, I want to have data reported for
#                                     # every position, whether it is present in the vcf file or not, so I will use count
#                                     # as the position, "-" as the seq, and 0 as the quality and depth of coverage
#                                     if int(pos) > count:
#                                         #print int(pos) - count, pos, count, range(count, int(pos))
#                                         # the number of skipped positions is equal to the value for pos - count
#                                         # For each skipped position (i), set appropriate variables to appropriate values
#                                         for i in range(count, int(pos)):
#                                             posAdj = count
#                                             seqAdj = "-"
#                                             matchAdj = 0
#                                             qualAdj = 0
#                                             DPAdj = 0
#                                             #vcfData[fileName][rL][fC][target][size][int(posAdj)][seqAdj][matchAdj][qualAdj] = DP
#                                             arrQual.append(float(qualAdj))
#                                             arrCov.append(float(DPAdj))
#                                             arrSum.append(float(matchAdj))
#                                             count += 1
#                                             if int(pos) == count:
#                                                 #vcfData[fileName][rL][fC][target][size][int(pos)][seq][match][qual] = DP
#                                                 arrQual.append(float(qual))
#                                                 arrCov.append(float(DP))
#                                                 arrSum.append(float(match))
#                                                 count += 1
#                                     else:
#                                         #vcfData[fileName][rL][fC][target][size][int(pos)][seq][match][qual] = DP
#                                         arrQual.append(float(qual))
#                                         arrCov.append(float(DP))
#                                         arrSum.append(float(match))
#                                         count += 1
#                         # In the case of no data being present in a file,
#                         total = count - 1
#                         if total == 0:
#                             avgQual = 0
#                             stdQual = 0
#                             avgCov = 0
#                             stdCov = 0
#                             avgID = 0
#                             qualMet = 0
#                         else:
#                             avgQual = sum(arrQual)/total
#                             stdQual = numpy.std(arrQual)
#                             avgCov = sum(arrCov)/total
#                             stdCov = numpy.std(arrCov)
#                             avgID = sum(arrSum)/total * 100
#                             qualMet = avgQual * avgCov
#
#                         outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
#                                     % (rLength, fCov, fileNoExt, size, avgQual, stdQual, avgCov, stdCov, avgID, qualMet))
#
#                         output.close()
#     outFile.close()


def run():
    global targetPath, sequencePath
    loadTargetData()
    targetFiles, sequenceFiles, sequenceNames = fileGrabber(targetPath, sequencePath)
    indexTargets(targetFiles)



# Run the pipeline
if __name__ == '__main__':
    run()
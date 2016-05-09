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

# Argument parser for user-inputted values, and a nifty help menu
from argparse import ArgumentParser

__author__ = 'adamkoziol'

# Parser for arguments
parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
parser.add_argument('-v', '--version', action='version', version='%(prog)s v1.0')
parser.add_argument('-p', '--path', required=True, help='Specify input directory')
parser.add_argument('-s', '--sequencepath', help='Path of .fastq(.gz) files to process. If not '
                    'provided, the default path of "path/sequences" will be used')
parser.add_argument('-t', '--targetpath', help='Path of target files to process. If not '
                    'provided, the default path of "path/targets" will be used')
parser.add_argument('-m', '--miSeqPath', help='Path of the folder containing MiSeq run data folder')
parser.add_argument('-f', '--miseqfolder', help='Name of the folder containing MiSeq run data')
parser.add_argument('-r1', '--readLengthForward', help='Length of forward reads to use. Can specify'
                    '"full" to take the full length of forward reads specified on the SampleSheet')
parser.add_argument('-r2', '--readLengthReverse', default=0, help='Length of reverse reads to use. '
                    'Can specify "full" to take the full length of reverse reads specified on the SampleSheet')
parser.add_argument('-c', '--customSampleSheet', help='Path of folder containing a custom sample '
                    'sheet (still must be named "SampleSheet.csv")')
parser.add_argument('-P', '--projectName', help='A name for the analyses. If nothing is provided, then '
                    'the "Sample_Project" field in the provided sample sheet will be used. Please note that bcl2fastq '
                    'creates subfolders using the project name, so if multiple names are provided, the results will be '
                    'split as into multiple projects')
parser.add_argument('-16S', '--16Styping', action='store_true', help='Perform 16S typing. Note that'
                    'for analyses such as MLST, pathotyping, serotyping, and virulence typing that require the genus'
                    'of a strain to proceed, 16S typing will still be performed')
parser.add_argument('-M', '--Mlst', action='store_true', help='Perform MLST analyses')
parser.add_argument('-Y', '--pathotYping', action='store_true', help='Perform pathotyping analyses')
parser.add_argument('-S', '--Serotyping', action='store_true', help='Perform serotyping analyses')
parser.add_argument('-V', '--Virulencetyping', action='store_true',
                    help='Perform virulence typing analyses')
parser.add_argument('-a', '--armi', action='store_true',
                    help='Perform ARMI antimicrobial typing analyses')
parser.add_argument('-r', '--rmlst', action='store_true', help='Perform rMLST analyses')
parser.add_argument('-d', '--detailedReports', action='store_true', help='Provide detailed reports with'
                    'percent identity and depth of coverage values rather than just "+" for positive results')
parser.add_argument('-C', '--customTargetPath', help='Provide the path for a folder of custom targets'
                    '.fasta format')

# TODO Add custom cutoffs
# TODO Assert .fastq files present in provided folder
# TODO Don't touch .fastq(.gz) files

# Get the arguments into a list
args = vars(parser.parse_args())

# Define variables from the arguments - there may be a more streamlined way to do this
path = os.path.join(args['path'], "")

# Optional path arguments. If these paths are present, add a trailing "/" with os.path.join
if args['miSeqPath']:
    miseqpath = os.path.join(args['miSeqPath'], "")
else:
    miseqpath = ""
if args['customSampleSheet']:
    customsamplesheet = os.path.join(args['customSampleSheet'], "")

# If these variables are not defined, set to default values. I chose to set the default values this way rather than
# specifying default=... in the parser.add_argument because I need to use variables that were defined following the
# argument parsing
if args['sequencepath']:
    sequencepath = os.path.join(args['sequencepath'], "")
else:
    sequencepath = "%ssequences/" % path

if args['targetpath']:
    targetpath = os.path.join(args['targetpath'], "")
else:
    targetpath = "%stargets/" % path

# Raw read sipping arguments
miseqfolder = args['miseqfolder']
projectname = args['projectName']
forwardreads = args['readLengthForward']
reversereads = args['readLengthReverse']


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    try:
        # os.makedirs makes parental folders as required
        os.makedirs(inpath)
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the dictionary responsible for storing the report data
metadatafile = defaultdict(make_dict)
# Define global variables
seqdict = defaultdict(make_dict)
# Create the name for the report folder - this variable will be passed on to other functions
reportfolder = "%sreports/%s" % (path, time.strftime("%Y.%m.%d.%H.%M.%S"))
make_path(reportfolder)
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


def printtime(string):
    """Prints a string in bold with the elapsed time
    :param string: a string to be printed in bold
    """
    global start
    print '\n\033[1m' + "[Elapsed Time: %0.2f seconds] %s" % (time.time() - start, string) + '\033[0m'


def foldererprepprocesses(samplename):
    """
    A helper function to make a pool of processes to allow for a multi-processed approach to error correction
    :param samplename: string of the name of the strain to decompress
    """
    global sequencepath
    print "Decompressing .gz files"
    # Initialise variables
    foldererprepargs = []
    createfoldererpool = Pool()
    # Prepare a tuple of the argument (name)
    for name in samplename:
        foldererprepargs.append(name)
    # This map function allows for multi-processing
    createfoldererpool.map(folderer, foldererprepargs)


def folderer(gzfile):
    """
    Uses gzip to decompress .gz.fastq files
    :param gzfile: string of file name and path of a gzipped file
    """
    # Define the gzip command line call
    gzipcommand = "gzip -d --force %s" % gzfile
    # Run the call
    subprocess.call(gzipcommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    dotter()


def filer(filelist):
    """
    Helper script that creates a set of the stain names created by stripping off parts of the filename.
    Hopefully handles different naming conventions (e.g. 2015-SEQ-001_S1_L001_R1_001.fastq(.gz),
    2015-SEQ-001_R1_001.fastq.gz, 2015-SEQ-001_R1.fastq.gz, 2015-SEQ-001_1.fastq.gz, and 2015-SEQ-001_1.fastq.gz
    all become 2015-SEQ-001)
    :param filelist: List of files to parse
    """
    # Initialise the set
    fileset = set()
    for seqfile in filelist:
        # Search for the conventional motifs present following strain names
        # _S\d+_L001_R\d_001.fastq(.gz) is a typical unprocessed Illumina fastq file
        if re.search("_S\d+_L001", seqfile):
            fileset.add(re.split("_S\d+_L001", seqfile)[0])
        # Files with _R\d_001.fastq(.gz) are created in the SPAdes assembly pipeline
        elif re.search("_R\d_001", seqfile):
            fileset.add(re.split("_R\d_001", seqfile)[0])
        # _R\d.fastq(.gz) represents a simple naming scheme for paired end reads
        elif re.search("R\d.fastq", seqfile):
            fileset.add(re.split("_R\d.fastq", seqfile)[0])
        # _\d.fastq is always possible
        elif re.search("[-_]\d.fastq", seqfile):
            fileset.add(re.split("[-_]\d.fastq", seqfile)[0])
        # .fastq is the last option
        else:
            fileset.add(re.split(".fastq", seqfile)[0])
        dotter()
    return fileset


def nameextractor():
    """Creates a set of the names of the files to be used in the analyses"""
    global sequencepath
    fileset = set()
    foldernames = set()
    filefull = []
    # Create lists of the .gz, the .fastq, and the folders in the path
    gzchecker = glob("%s*.gz" % sequencepath)
    # Extracts the .gz files as required
    if gzchecker:
        foldererprepprocesses(gzchecker)
    fastqchecker = glob("%s*.fastq" % sequencepath)
    # For each list, ensure that the list exists...
    if fastqchecker:
        # Create appropriately named folders and move the .fastq files into the proper folder
        filelist, foldernames, filefull = seqmovr(fastqchecker)
    folderchecker = glob("%s*/" % sequencepath)
    # Get folder names in the sequences folder - this should only be useful if there were some files
    # in folders, and some not in folders
    if folderchecker:
        for seqname in folderchecker:
            filefull = []
            # Add the appropriate variables to the appropriate sets
            fileset.add(seqname)
            foldernames.add(os.path.split(seqname)[0].split("/")[-1])
            for fastqfile in glob("%s/*.fastq" % seqname):
                # Baited files shouldn't be added to the list now.
                if not re.search("bait_match", fastqfile):
                    filefull.append(fastqfile)
            dotter()
            #  Add the unique fastq files to the dictionary
            fastqfnr = []
            foldername = os.path.split(seqname)[0].split("/")[-1]
            for files in filefull:
                # Use fileR to determine the "short (and unique) file name"
                filename = filer([os.path.basename(files)]).pop()
                # If these variables are identical, add the file to the list
                if foldername == filename:
                    fastqfnr.append(files)
                # Add the list to the dictionary
                seqdict[foldername]["rawFastq"] = sorted(fastqfnr)
    # Return the lists of path/folder names and folder names
    return list(fileset), list(foldernames), filefull


def seqmovr(seqfull):
    """
    Creates appropriately named folders, and moves sequence files to appropriate folders
    :param seqfull: list of all strains in the analysis
    """
    global sequencepath, count
    count = 0
    foldernames = set()
    filefull = []
    # Extract the unique names using the helper function fileR...
    filelist = filer(seqfull)
    # Create the folders if necessary
    print "\nMoving files to appropriate folders"
    for folder in filelist:
        # There can be an empty string in this set (''), which breaks the loop
        if folder:
            # Split the path from the folders - will give a name like 2015-SEQ-794
            foldername = os.path.split(folder)[-1]
            foldernames.add(foldername)
            #  Make folders as necessary
            make_path(folder)
            # Search the path for any file or folder that contains the seqName
            filecheck = [f for f in os.listdir(sequencepath) if re.search(foldername, f)]
            for seqfile in filecheck:
                # Use fileR to determine the "short" name of the sequence file
                seqset = filer([seqfile]).pop()
                # Move files (ignore folders) if folderName identical to the "short" file name
                if os.path.isfile("%s%s" % (sequencepath, seqfile)) and foldername == seqset:
                    shutil.move("%s%s" % (sequencepath, seqfile), "%s%s/%s" % (sequencepath, foldername, seqfile))
                    filefull.append("%s%s/%s" % (sequencepath, foldername, seqfile))
    return filelist, foldernames, filefull


def baitrprocesses(analysistype):
    """
    Baiting multiprocessing helper function
    :param analysistype: string of the analysis type
    """
    global seqdict
    # Initialise the args list
    baitargs = []
    # Initialise the pool of processes - it defaults to the number of processors
    baitpool = Pool()
    for foldername in seqdict:
        baittype = seqdict[foldername]["bait"]["fastqFiles"]
        # Iterate through the different bait types stored in the dictionary
        for btype in baittype.keys():
            #  Only perform the baiting for the current analysis type
            if analysistype in btype:
                # Append the appropriate variables
                baitargs.append((seqdict[foldername]["bait"]["fastqFiles"][btype],
                                 seqdict[foldername]["rawFastq"], btype))
    # Map the baitArgs to the baitR function
    folderpath = (baitpool.map(baitr, baitargs))
    # Repopulate the sequence dictionary with the new bait files
    for folder in folderpath:
        foldername = os.path.split(os.path.split(folder)[0])[1]
        baitmatch = glob("%s/*.fastq" % folder)
        baittype = seqdict[foldername]["bait"]["fastqFiles"].keys()[0]
        seqdict[foldername]["bait"]["fastqFiles"][baittype] = baitmatch


def baitr((baitfile, sequences, baittype)):
    """Runs mirabait with selected targets on fastq files"""
    global targetpath
    # Works on paired-end and unpaired-end reads
    if len(sequences) == 2:
        # Define the forward and reverse .fastq reads based on their position in the (ordered) list
        forwardfastqpath = sequences[0]
        reversefastqpath = sequences[1]
        # Get the path of the fastq files
        baitpath = os.path.split(forwardfastqpath)[0] + "/%s" % baittype
        # Create the path (if necessary)
        make_path(baitpath)
        baitcheck = glob("%s/%s_match*" % (baitpath, baittype))
        # Check to see if the baiting process has been performed previously
        if len(baitcheck) < 2:
            # Define the baiting command - if a hash of the target has previously been computed and placed in the
            # bait folder, then proceed appropriately
            if ".gz" in baitfile:
                baitcommand = "cd %s && mirabait -L %s -N %s -p %s %s" \
                              % (baitpath, baitfile, baittype, forwardfastqpath, reversefastqpath)
                # Define /dev/null
                fnull = open(os.devnull, 'wb')
                # Run the command
                subprocess.call(baitcommand, shell=True, stdout=fnull, stderr=fnull)
            else:
                baitcommand = "cd %s && mirabait -b %s -N %s -p %s %s" \
                              % (baitpath, baitfile, baittype, forwardfastqpath, reversefastqpath)
                # Move the newly created hashstat file to the bait folder
                baithash = "%s/hashstat.mhs.gz" % baitpath
                hashdestination = baitfile.split(".")[0] + ".mhs.gz"
                # Define /dev/null
                fnull = open(os.devnull, 'wb')
                # Run the command
                subprocess.call(baitcommand, shell=True, stdout=fnull, stderr=fnull)
                # Try to move the baitHash file to the bait path, but pass on failure
                try:
                    shutil.move(baithash, hashdestination)
                except IOError:
                    pass
    else:
        # Define the forward and reverse .fastq reads based on their position in the (ordered) list
        forwardfastqpath = sequences[0]
        # Get the path of the fastq files
        baitpath = os.path.split(forwardfastqpath)[0] + "/%s" % baittype
        make_path(baitpath)
        baitcheck = glob("%s/%s_match*" % (baitpath, baittype))
        if not baitcheck:
            # Define the baiting command
            if ".gz" in baitfile:
                baitcommand = "cd %s && mirabait -L %s -N %s %s" % (baitpath, baitfile, baittype, forwardfastqpath)
                # Define /dev/null
                fnull = open(os.devnull, 'wb')
                # Run the command
                subprocess.call(baitcommand, shell=True, stdout=fnull, stderr=fnull)
            else:
                baitcommand = "cd %s && mirabait -b %s -N %s %s" % (baitpath, baitfile, baittype, forwardfastqpath)
                # Define /dev/null
                fnull = open(os.devnull, 'wb')
                # Run the command
                subprocess.call(baitcommand, shell=True, stdout=fnull, stderr=fnull)
                # Move the newly created hashstat file to the bait folder
                # Set the hash name and the destination
                baithash = "%s/hashstat.mhs.gz" % baitpath
                hashdestination = baitfile.split(".")[0] + ".mhs.gz"
                # Try to move the file, but pass on an IOError
                try:
                    shutil.move(baithash, hashdestination)
                except IOError:
                    pass
    dotter()
    return baitpath


def sixteensreportmaker(resultdict, analysistype):
    """
    Creates reports for 16S analyses
    :param resultdict: dictionary containing results of 16S analyses
    :param analysistype: string of the analysis type
    """
    global reportfolder, seqdict
    # Create the path if necessary
    make_path(reportfolder)
    # Initialise variables
    compiledresultstring = ""
    csvheader = ""
    # Iterate through each strain in resultDict
    for strain in seqdict:
        # Create variables as required
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        reportdir = "%s/reports" % fastqdir
        # Create a strain-specific report in the <strain>/16S/reports directory
        reportname = "%s/%s_%s_reports.tsv" % (reportdir, strain, baittype)
        make_path(reportdir)
        csvfile = open(reportname, "wb")
        # Get the header started
        csvheader = "Strain,Genus,Accession,PercentIdentity,AvgFoldCov\n"
        # Initialise the variable to hold the result data
        resultstring = ""
        # The results are stored in a nested dictionary
        for genus in resultdict[strain]:
            for database in resultdict[strain][genus]:
                for accession in resultdict[strain][genus][database]:
                    for avgidentity, avgfoldcov in resultdict[strain][genus][database][accession].iteritems():
                        # Add the appropriate variables to the string
                        resultstring = "%s,%s,%s,%s,%s" % (strain, genus, accession, avgidentity, avgfoldcov)
        # Append the result plus a newline character to string containing results for all strains
        compiledresultstring += resultstring + "\n"
        # Write the header and data strings to file
        csvfile.write(csvheader)
        csvfile.write(resultstring)
        csvfile.close()
    # Create a report containing all the 16S results from each strain
    compiledcsvfile = open("%s/%s_results.tsv" % (reportfolder, analysistype), "wb")
    # Write the header and data strings to file
    compiledcsvfile.write(csvheader)
    compiledcsvfile.write(compiledresultstring)
    compiledcsvfile.close()


def pathoreportr(matchdict, analysistype, organismdict, organismlist):
    """
    Creates reports for pathotyping and serotyping results
    :param matchdict: dictionary of hits of a strain to targets
    :param analysistype: string of the analysis type
    :param organismdict: dictionary of the 16S results
    :param organismlist: list of all the genera in the current analysis
    """
    global reportfolder, seqdict, args
    # Initialise data string
    completestring = ""
    # Iterate through the genera present in the analyses
    for organism in organismlist:
        # Initialise variables
        orgcount = 0
        organismdata = ""
        organismheader = "Strain,Genus,"
        # Iterate through all the strains in the analysis
        for strain in seqdict:
            # Try/except for attribute errors - if a strain does not have results, then just pass
            try:
                strainheader = "Strain,"
                straindata = ""
                # Iterate through the organisms in organismdict
                for currentorganism in organismdict[strain]:
                    # If this current organism matches the organism of interest, proceed
                    if organism == currentorganism:
                        # Increment orgcount
                        orgcount += 1
                        headergenes = ''
                        # Iterate through the targets to find the names of the genes for the headers
                        for target in sorted(seqdict[strain]["targets"][analysistype]):
                            # Set target name
                            targetname = os.path.basename(target).split(".")[0]
                            headergenes += "%s," % targetname
                        # As the header and data variables are strings, the header information only needs to be appended
                        # once, otherwise, the header would precede the data for every strain
                        if orgcount == 1:
                            organismheader += headergenes
                            completestring += "Strain,Genus,"
                            completestring += headergenes
                        # Strain header is unique for each strain, so it should be populated every time
                        strainheader += headergenes
                        # Populate the data strings with the appropriate variables
                        organismdata += "\n%s,%s," % (strain, organism)
                        completestring += "\n%s,%s," % (strain, organism)
                        # Organism is not necessary for the strain-specific reports
                        straindata += "\n%s," % strain
                        # Iterate through the targets to find the presence/absence of each target
                        for target in sorted(seqdict[strain]["targets"][analysistype]):
                            targetname = os.path.basename(target).split(".")[0]
                            # If there is a result for a strain/target combination
                            if matchdict[strain][targetname]:
                                # Pathotype results are treated slightly differently than serotyping results
                                if analysistype == "pathotype":
                                    # Allow for the inclusion of more data in the results table
                                    if args['detailedReports']:
                                        gene = matchdict[strain][targetname].keys()[0]
                                        percentidentity = matchdict[strain][targetname][gene].keys()[0]
                                        foldcoverage = matchdict[strain][targetname][gene].values()[0]
                                        completestring += "%s;%s;%s," % (gene, percentidentity, foldcoverage)
                                        organismdata += "%s;%s;%s," % (gene, percentidentity, foldcoverage)
                                        straindata += "%s;%s;%s," % (gene, percentidentity, foldcoverage)
                                    else:
                                        completestring += "+,"
                                        organismdata += "+,"
                                        straindata += "+,"

                                else:
                                    # Serotyping is not a binary presence/absence. The calculated serotype is reported
                                    serotype = matchdict[strain][targetname].keys()[0].split("_")[-1]
                                    # Populate the strings with the serotype
                                    completestring += "%s," % serotype
                                    organismdata += "%s," % serotype
                                    straindata += "%s," % serotype
                            # If there is no strain/target result, populate the strings with negative results ("-")
                            else:
                                completestring += "-,"
                                organismdata += "-,"
                                straindata += "-,"
                    # Pull values from seqdict in order to create properly named reports in the appropriate folders
                    baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
                    fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
                    reportdir = "%s/reports" % fastqdir
                    reportname = "%s/%s_%s_reports.csv" % (reportdir, strain, baittype)
                    make_path(reportdir)
                    # Open the strain-specific report, and write the appropriate header/results
                    straincsvfile = open(reportname, "wb")
                    straincsvfile.write(strainheader)
                    straincsvfile.write(straindata)
                    straincsvfile.close()
                # Create the organism-specific report, and write the appropriate header/results
                organismcsvfile = open("%s/%s_%s_results.csv" % (reportfolder, analysistype, organism), "wb")
                organismcsvfile.write(organismheader)
                organismcsvfile.write(organismdata)
                organismcsvfile.close()
            # If there are no pathotyping/serotyping results for a particular strain, pass
            except AttributeError:
                pass
        # Add a newline character to completestring for each strain
        completestring += "\n"
    # Create the report for each strain in the data set, and write the appropriate header/results
    compiledcsvfile = open("%s/%s_results.csv" % (reportfolder, analysistype), "wb")
    compiledcsvfile.write(completestring)
    compiledcsvfile.close()


def virulencer(resultsdict, analysistype):
    """
    Determines genus of strains based on reference mapping results
    :param resultsdict: dictionary containing the best match to a target
    :param analysistype: string of analysis type
    """
    global seqdict
    notefile = ""
    strainheader = ""
    combinedvirulence = []
    detailedresults = ""
    # Iterate through the strains in resultsdict
    for strain in seqdict:
        # Initialise the set containing the virulence types
        virulenceset = set()
        # Retrieve the target path from seqdict
        for target in seqdict[strain]["targets"][analysistype]:
            targetspath = os.path.split(target)[0]
            # The notes.txt file stores the associations between gene names and subtype
            notefile = open(targetspath + "/notes.txt", "rb").readlines()
        # Iterate through the nested dictionary
        for targetname in resultsdict[strain]:
            for allele in resultsdict[strain][targetname]:
                for line in notefile:
                    # The format of the gene name is something like stx2:3:GQ429162:3, while the corresponding entry
                    # in notes.txt would be: stx23:ONT 23765, variant a.  Therefore targetname (stx2) +
                    # allele.split(":")[1] (3) + ":" yields stx23: If that is in the line, proceed
                    alleleresult = allele.split(":")[0] + allele.split(":")[-1]
                    if alleleresult in line:
                        # if targetname + allele.split(":")[1] + ":" in line:
                        # Add the target name (stx2) + the trailing two characters in the line (d1) to the set
                        if resultsdict[strain][targetname][allele].keys()[0] >= 98.0:
                            # print strain, alleleresult, alleleresult[4:5], resultsdict[strain][targetname][allele].keys()[0]
                            virulenceset.add(alleleresult)
                        # There are a two stx2 variants with a subtype: d1 and d2. Allow for this exception
                        # if line.rstrip()[-2] == "d":
                        #
                        #     # Create the string for the detailed reports
                        #     detailedresults += "%s,%s,%s,%s\n" % (strain, alleleresult + line.rstrip()[-2],
                        #                                           resultsdict[strain][targetname][allele].keys()[0],
                        #                                           resultsdict[strain][targetname][allele].values()[0])
                        # else:
                        # Add the target name (stx2) + the trailing character (a) in the line -> (stx2a) to the set
                        detailedresults += "%s,%s,%s,%s\n" % (strain, alleleresult,
                                                              resultsdict[strain][targetname][allele].keys()[0],
                                                              resultsdict[strain][targetname][allele].values()[0])
        # Add an extra line to the detailed results to delineate between strains
        detailedresults += "\n"
        # Pull values from seqdict in order to create properly named reports in the appropriate folders
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        reportdir = "%s/reports" % fastqdir
        reportname = "%s/%s_%s_reports.csv" % (reportdir, strain, baittype)
        make_path(reportdir)
        # Open the strain-specific report, and write the appropriate header/results
        strainheader = "Strain,virulenceType\n"
        vtxset = set()
        for vtx in sorted(list(virulenceset)):
            vtype = vtx[:4] + vtx[-1]
            # subunit = vtx[4:5]
            if vtype != 'stx2b':
                if vtx[:4] + 'A' + vtx[-1] in virulenceset and vtx[:4] + 'B' + vtx[-1] in virulenceset:
                    # print strain, vtx, vtype, subunit, vtx[:4] + 'A' + vtx[-1]
                    vtxset.add(vtype)
            else:
                # print strain, '2b'
                vtxset.add(vtype)

        straindata = "%s,%s\n" % (strain, "; ".join(sorted(list(vtxset))))
        combinedvirulence.append(straindata)
        straincsvfile = open(reportname, "wb")
        straincsvfile.write(strainheader)
        straincsvfile.write(straindata)
        straincsvfile.close()
    # Create the organism-specific report, and write the appropriate header/results
    organismcsvfile = open("%s/%s_results.csv" % (reportfolder, analysistype), "wb")
    organismcsvfile.write(strainheader)
    organismcsvfile.write("".join(combinedvirulence))
    organismcsvfile.close()
    if args['detailedReports']:
        detailedheader = "Strain,virulenceType,PercentIdentity,AverageCoverage\n"
        detailedreport = open("%s/%s_DetailedResults.csv" % (reportfolder, analysistype), "wb")
        detailedreport.write(detailedheader)
        detailedreport.write(detailedresults)
        detailedreport.close()


def reportr(inputdict, analysistype):
    """
    Creates reports for custom analyses (and perhaps future, generic analyses as well?)
    :param inputdict: dictionary containing genesipping results to be formatted into a report
    :param analysistype: string of the name of the current analysis
    """
    global seqdict, reportfolder
    header = ""
    sippingresults = ""
    detailedresults = ""
    # Iterate through the strain with results in inputdict
    for strain in seqdict:
        strainresults = "%s," % strain
        sippingresults += "%s," % strain
        detailedresults += "%s," % strain
        straindetails = "%s," % strain
        # Populate the header
        header = "Strain,"
        for target in sorted(seqdict[strain]["targets"][analysistype]):
            targetname = os.path.basename(target).split(".")[0]
            header += targetname + ","
        # Iterate through the targets in seqdict
        for target in sorted(seqdict[strain]["targets"][analysistype]):
            targetname = os.path.basename(target).split(".")[0]
            # If there are results in the dictionary, populate the appropriate variables
            if inputdict[strain][targetname]:
                for results in inputdict[strain][targetname]:
                    # Support for multifasta files - if the name of the gene doesn't match the file name, print the
                    # gene name instead of a +
                    if results == targetname:
                        sippingresults += "+,"
                        strainresults += "+,"
                        detailedresults += "%s;%s," % (inputdict[strain][targetname][results].keys()[0],
                                                       inputdict[strain][targetname][results].values()[0])
                        straindetails += "%s;%s," % (inputdict[strain][targetname][results].keys()[0],
                                                     inputdict[strain][targetname][results].values()[0])
                        # Catches multiple hits - only keeps the best result
                        break
                    else:
                        sippingresults += "%s," % results
                        strainresults += "%s," % results
                        detailedresults += "%s;%s;%s," % (results, inputdict[strain][targetname][results].keys()[0],
                                                          inputdict[strain][targetname][results].values()[0])
                        straindetails += "%s;%s;%s," % (results, inputdict[strain][targetname][results].keys()[0],
                                                        inputdict[strain][targetname][results].values()[0])
                        # Catches multiple hits - only keeps the best result
                        break
            # Populate the string with - to indicate the target was not found
            else:
                sippingresults += "-,"
                strainresults += "-,"
                detailedresults += "-,"
                straindetails += "-,"
        # Add newline characters for proper formatting
        sippingresults += "\n"
        detailedresults += "\n"
        header += "\n"
        # Pull values from seqdict in order to create properly named reports in the appropriate folders
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        reportdir = "%s/reports" % fastqdir
        reportname = "%s/%s_%s_reports.csv" % (reportdir, strain, baittype)
        make_path(reportdir)
        with open(reportname, "wb") as report:
            report.write(header)
            if args['detailedReports']:
                report.write(straindetails)
            else:
                report.write(strainresults)
    # Create the report for each strain in the data set, and write the appropriate header/results
    with open("%s/%s_results.csv" % (reportfolder, analysistype), "wb") as report:
        report.write(header)
        if args['detailedReports']:
            report.write(detailedresults)
        else:
            report.write(sippingresults)


def baittargets(currenttargetpath, analysistype):
    """
    Creates a file to be used for baiting if one does not already exist by concatenating together all individual
    target files in the current target path
    :param currenttargetpath: path of the targets in the current analysis
    :param analysistype: string of the current analysis type
    """
    if not os.path.isfile("%s/bait/%sBait.fa" % (currenttargetpath, analysistype)):
        # Creates the bait path if necessary
        make_path("%s/bait" % currenttargetpath)
        # Get a list of the targets into a variable - use a list comprehension to remove any database files present
        # target.fa, or target.fasta should be acceptable, but target.fa.nsi, or target.fasta.smi should be eliminated
        targetfiles = [fasta for fasta in glob("%s/*.fa*" % currenttargetpath) if not re.search("\.fa.*[\.]", fasta) and
                       "Concatenated" not in fasta]
        # Initialise a list to store all the target records
        alltargets = []
        # Iterate through the filtered list
        for tfile in targetfiles:
            # Pull all records using BioPython
            for record in SeqIO.parse(open(tfile, "rU"), "fasta"):
                # Append the records to the list
                alltargets.append(record)
        # Create the name of the output file
        concatenated = "%s/bait/%sBait.fa" % (currenttargetpath, analysistype)
        # Open the file
        with open(concatenated, "wb") as formatted:
            # Write all the records to the concatenated file
            SeqIO.write(alltargets, formatted, "fasta")


def sixteens():
    """Performs the necessary analyses on strains using 16S targets"""
    global targetpath, seqdict
    # Set the analysis type variable to 16S. This variable is important for retrieving 16S-specific data from seqdict
    analysistype = "16S"
    # Set the path of the analysistype data
    currenttargetpath = "%s%s" % (targetpath, analysistype)
    # Create the bait target file (if necessary)
    baittargets(currenttargetpath, analysistype)
    # In order to save time, during baiting, a precomputed hash file is used
    hashfile = glob("%s/bait/*.gz" % currenttargetpath)
    # If the hashfile is present, use it
    if hashfile:
        baitfile = hashfile[0]
    # Otherwise, use the concatenated bait file
    else:
        baitfile = glob("%s/bait/*.fa*" % currenttargetpath)[0]
    print "Filtering .fastq files with %s targets" % analysistype
    # Get a list of the targets into a variable - use a list comprehension to remove any database files present
    # target.fa, or target.fasta should be acceptable, but target.fa.nsi, or target.fasta.smi should be eliminated
    sixteensdb = [fasta for fasta in glob("%s/*.fa*" % currenttargetpath) if not re.search("\.fa.*[\.]", fasta) and
                  "Concatenated" not in fasta]
    # Populate seqdict with the name and path of the bait file used as well as the database
    for foldername in seqdict:
        seqdict[foldername]["bait"]["fastqFiles"][analysistype] = baitfile
        seqdict[foldername]["targets"][analysistype] = sixteensdb
    # Run the baiting process
    baitrprocesses(analysistype)
    # Perform SMALT indexing of targets
    print "\nIndexing %s targets" % analysistype
    SMALT.smaltindextargetsprocesses(sixteensdb, currenttargetpath)
    # Perform reference mapping with SMALT
    print '\nPerforming %s reference mapping' % analysistype
    SMALT.smaltmappingprocesses(seqdict, analysistype, "SMALT")
    # Use samtools to sort the bam files
    print "\nSorting mapped %s files" % analysistype
    bamProcessor.sortingprocesses(seqdict, analysistype)
    # Use samtools to index the sorted bam files
    print '\nIndexing sorted %s files' % analysistype
    bamProcessor.bamindexingprocesses(seqdict, analysistype)
    # Use pysamstats to parse the sorted, indexed bam files
    print '\nParsing %s results' % analysistype
    generadict, generalist = bamPysamStats.bamparseprocesses(seqdict, analysistype)
    # Create reports of the results
    sixteensreportmaker(generadict, analysistype)
    # Return the computed genera
    return generadict, generalist


def mlst(organismdict, organismlist):
    """
    Performs the necessary analyses on strains using genus-specific MLST targets
    :param organismdict: dictionary of the 16S results
    :param organismlist: list of all the genera in the current analysis
    """
    global targetpath, seqdict
    # Set the analysis type
    analysistype = "MLST"
    profiletype = "MLSTprofile"
    # MLST targets are stored in targetpath/Organism/<genus>/MLST/alleles
    currenttargetpath = "%sOrganism" % targetpath
    # Print this out before the loop
    print 'Indexing %s target files' % analysistype
    for strain in organismdict:
        # If there is not an MLST scheme installed for a particular organism, then the script will crash when it tries
        # to find the necessary files, as they are not present. Allow index errors to pass
        try:
            # Using the organismdict entry (genus) generated in the 16S analysis, set the allele and profile paths
            # NB: This will come up multiple times with this script, but I only allowed a small amount of freedom in
            # the placement of folders. Usually, there are strict folder hierarchies, which must be followed
            mlstpath = glob("%s/%s/*MLST*" % (currenttargetpath, organismdict[strain].keys()[0]))[0] + "/alleles"
            profilepath = glob("%s/%s/*MLST*" % (currenttargetpath, organismdict[strain].keys()[0]))[0] + "/profile"
            # Try to find the precomputed .json profile file
            profilefile = glob("%s/*.json" % profilepath)
            # If it does not exist, find the .txt profile file. As the script requires a small amount of formatting on
            # the profile file prior to analysis, I forced a changed in file extension to hopefully ensure that this
            # formatting has been performed
            if not profilefile:
                profilefile = glob("%s/*.txt" % profilepath)[0]
            else:
                profilefile = profilefile[0]
            # Create the bait target file (if necessary)
            baittargets(currenttargetpath, analysistype)
            # If a precomputed hash file is present in the bait folder, use it to save on processing time
            hashfile = glob("%s/bait/*.gz" % mlstpath)
            if hashfile:
                baitfile = hashfile[0]
            # Otherwise, use the .fasta bait file created above
            else:
                baitfile = glob("%s/bait/*.fa*" % mlstpath)[0]
            # Set the bait type variable using the genus of the strain and the analysis type
            baittype = "%s_MLST" % organismdict[strain].keys()[0]
            # Store the baittype variable in seqdict
            seqdict[strain]["bait"]["fastqFiles"][baittype] = baitfile
            # Get a list of the targets into a variable - use a list comprehension to remove any database files present
            # target.fa, or target.fasta should be acceptable, but target.fa.nsi, or target.fasta.smi will be eliminated
            targets = [fasta for fasta in glob("%s/*.fa*" % mlstpath) if not re.search("\.fa.*[\.]", fasta) and
                       "Concatenated" not in fasta]
            # Enter the targets and profile file into the dictionary
            seqdict[strain]["targets"][analysistype] = targets
            seqdict[strain]["MLSTprofile"] = profilefile
            # Index the SMALT targets
            SMALT.smaltindextargetsprocesses(targets, mlstpath)
        except IndexError:
            pass
    # Bait!
    print "\nFiltering .fastq files with %s targets" % analysistype
    baitrprocesses(analysistype)
    # Get the MLST profile into a dictionary
    print "\nLoading %s profiles" % analysistype
    profiledict, genedict = rawMLST.profilR(seqdict, profiletype)
    print '\nPerforming %s reference mapping' % analysistype
    # Perform SMALT reference mapping
    SMALT.smaltmappingprocesses(seqdict, analysistype, "SMALT")
    print "\nSorting mapped %s files" % analysistype
    # Use samtools to sort reference mapped bam files
    bamProcessor.sortingprocesses(seqdict, analysistype)
    print '\nIndexing sorted %s files' % analysistype
    # Use samtools to index sorted reference mapped bam files
    bamProcessor.bamindexingprocesses(seqdict, analysistype)
    print '\nParsing %s results' % analysistype
    # Use pysamstats to parse indexed sorted reference mapped bam files
    mlstmatches = bamPysamStats.bamparseprocesses(seqdict, analysistype)
    print '\nFinding multilocus sequence types'
    # Determine sequence types
    sequencetypes = rawMLST.sequenceTyper(mlstmatches, profiledict, genedict)
    # Create a report
    rawMLST.MLSTreportMaker(seqdict, sequencetypes, analysistype, reportfolder, organismdict, organismlist, path)


def pathotyper(organismdict, organismlist, analysistype):
    """
    Performs the necessary analyses on strains using genus-specific pathotype/serotype targets
    :param organismdict: dictionary of the 16S results
    :param organismlist: list of all the genera in the current analysis
    :param analysistype: string of the current analysis type
    """
    global targetpath, seqdict
    # Targets are stored in targetpath/Organism/<genus>/<analysistype>/
    currenttargetpath = "%sOrganism" % targetpath
    # Print this prior to the loop
    print 'Indexing %s target files' % analysistype
    # As strains will be processed depending differently based on their genus, they must be processed separately
    for strain in organismdict:
        # Try/except account for genera without pathotyping schemes
        try:
            # Set the target path
            pathopath = glob("%s/%s/%s" % (currenttargetpath, organismdict[strain].keys()[0], analysistype))[0]
            # Create the bait targets if necessary
            baittargets(pathopath, analysistype)
            # The cat file will be used in the "combined" reference mapping
            catfile = "%s/%sConcatenated.fasta" % (pathopath, analysistype)
            # If the cat file does not exist, copy the bait file from the bait folder to the target folder
            if not os.path.isfile(catfile):
                shutil.copyfile("%s/bait/%sBait.fa" % (pathopath, analysistype), catfile)
            # Find the files to be used in baiting - if a precomputed hash file is present, use it
            hashfile = glob("%s/bait/*.gz" % pathopath)
            if hashfile:
                baitfile = hashfile[0]
            else:
                baitfile = glob("%s/bait/*.fa*" % pathopath)[0]
            # Set the baittype variable as the genus, and the analysis type
            baittype = "%s_%s" % (organismdict[strain].keys()[0], analysistype)
            # Get all the targets into a list
            # Even though the cat file will be used for the reference mapping rather than the individual target files,
            # the target names (taken from the name of the target files) is still used in the parsing of results
            # Use a list comprehension to remove any database files present: target.fa, or target.fasta should be
            # acceptable, but target.fa.nsi, or target.fasta.smi should be eliminated
            targets = [fasta for fasta in glob("%s/*.fa*" % pathopath) if not re.search("\.fa.*[\.]", fasta) and
                       "Concatenated" not in fasta]
            # Add the bait file, the cat file, and the list of targets to seqdict
            seqdict[strain]["bait"]["fastqFiles"][baittype] = baitfile
            seqdict[strain]["targets"][analysistype] = targets
            seqdict[strain]["concatenatedTargets"][analysistype] = catfile
            # Index the SMALT targets
            SMALTcombined.smaltindextargets(catfile, pathopath)
        except IndexError:
            pass
    # Bait!
    print "\nFiltering .fastq files with %s targets" % analysistype
    baitrprocesses(analysistype)
    print '\nPerforming reference mapping'
    # Use SMALT to perform reference mapping
    SMALTcombined.smaltmappingprocesses(seqdict, analysistype, "SMALT")
    print '\nSorting mapped %s files' % analysistype
    # Use samtools to sort bam files
    bamProcessorCombined.sortingprocesses(seqdict, analysistype)
    print '\nIndexing sorted %s files' % analysistype
    # Use samtools to index sorted bam files
    bamProcessorCombined.bamindexingprocesses(seqdict, analysistype)
    print '\nParsing %s results' % analysistype
    # Use pysamstats to parse results
    pathomatches = bamPysamStatsCombined.bamparseprocesses(seqdict, analysistype)
    # Create a report
    pathoreportr(pathomatches, analysistype, organismdict, organismlist)


def virulencetyper(organismdict, analysistype):
    """
    Performs the necessary analyses on strains using genus-specific pathotype/serotype targets
    :param organismdict: dictionary of the 16S results
    :param analysistype: string of the current analysis type
    """
    global targetpath, seqdict
    # Targets are stored in targetpath/Organism/<genus>/<analysistype>/
    currenttargetpath = "%sOrganism" % targetpath
    # Print this prior to the loop
    print 'Indexing %s target files' % analysistype
    # As strains will be processed depending differently based on their genus, they must be processed separately
    for strain in organismdict:
        # Try/except account for genera without pathotyping schemes
        try:
            # Set the target path
            viropath = glob("%s/%s/%s" % (currenttargetpath, organismdict[strain].keys()[0], analysistype))[0]
            # Create the bait targets if necessary
            baittargets(viropath, analysistype)
            # The cat file will be used in the "combined" reference mapping
            catfile = "%s/%sConcatenated.fasta" % (viropath, analysistype)
            # If the cat file does not exist, copy the bait file from the bait folder to the target folder
            if not os.path.isfile(catfile):
                shutil.copyfile("%s/bait/%sBait.fa" % (viropath, analysistype), catfile)
            # Find the files to be used in baiting - if a precomputed hash file is present, use it
            hashfile = glob("%s/bait/*.gz" % viropath)
            if hashfile:
                baitfile = hashfile[0]
            else:
                baitfile = glob("%s/bait/*.fa*" % viropath)[0]
            # Set the baittype variable as the genus, and the analysis type
            baittype = "%s_%s" % (organismdict[strain].keys()[0], analysistype)
            # Get all the targets into a list
            # Even though the cat file will be used for the reference mapping rather than the individual target files,
            # the target names (taken from the name of the target files) is still used in the parsing of results
            # Use a list comprehension to remove any database files present: target.fa, or target.fasta should be
            # acceptable, but target.fa.nsi, or target.fasta.smi should be eliminated
            targets = [fasta for fasta in glob("%s/*.fa*" % viropath) if not re.search("\.fa.*[\.]", fasta) and
                       "Concatenated" not in fasta]
            # Add the bait file, the cat file, and the list of targets to seqdict
            seqdict[strain]["bait"]["fastqFiles"][baittype] = baitfile
            seqdict[strain]["targets"][analysistype] = targets
            seqdict[strain]["concatenatedTargets"][analysistype] = catfile
            # Index the SMALT targets
            SMALT.smaltindextargetsprocesses(targets, viropath)
        except IndexError:
            pass
    # Bait!
    print "\nFiltering .fastq files with %s targets" % analysistype
    baitrprocesses(analysistype)
    print '\nPerforming reference mapping'
    # Use SMALT to perform reference mapping
    SMALT.smaltmappingprocesses(seqdict, analysistype, "SMALT")
    print '\nSorting mapped %s files' % analysistype
    # Use samtools to sort bam files
    bamProcessor.sortingprocesses(seqdict, analysistype)
    print '\nIndexing sorted %s files' % analysistype
    # Use samtools to index sorted bam files
    bamProcessor.bamindexingprocesses(seqdict, analysistype)
    print '\nParsing %s results' % analysistype
    # Use pysamstats to parse results
    viromatches = bamPysamStats.bamparseprocesses(seqdict, analysistype)
    # Create a report
    virulencer(viromatches, analysistype)


def armi():
    """Performs the necessary analyses on strains using armi targets"""
    global targetpath, seqdict
    # Set the analysis type
    analysistype = "ARMI"
    # Set the path of the target files
    currenttargetpath = "%s%s" % (targetpath, analysistype)
    # Create the bait files as necessary
    baittargets(currenttargetpath, analysistype)
    # The cat file is all the target files concatenated together
    catfile = "%s/%sConcatenated.fasta" % (currenttargetpath, analysistype)
    # If the cat file doesn't exist, create it from the bait file generated in baittargets()
    if not os.path.isfile(catfile):
        shutil.copyfile("%s/bait/%sBait.fa" % (currenttargetpath, analysistype), catfile)
    # In order to save time, a precomputed hash file is used
    hashfile = glob("%s/bait/*.gz" % currenttargetpath)
    # If this precomputed hash exists, use it
    if hashfile:
        baitfile = hashfile[0]
    # Otherwise use the .fasta bait file
    else:
        baitfile = glob("%s/bait/*.fa*" % currenttargetpath)[0]
    # Get a list of the targets into a variable - use a list comprehension to remove any database files present
    # target.fa, or target.fasta should be acceptable, but target.fa.nsi, or target.fasta.smi should be eliminated
    armidatabase = [fasta for fasta in glob("%s/*.fa*" % currenttargetpath) if not re.search("\.fa.*[\.]", fasta) and
                    "Concatenated" not in fasta]
    # Add necessary variables to seqdict
    for folderName in seqdict:
        seqdict[folderName]["bait"]["fastqFiles"][analysistype] = baitfile
        seqdict[folderName]["targets"][analysistype] = armidatabase
        seqdict[folderName]["concatenatedTargets"][analysistype] = catfile
    print "Filtering .fastq files with %s targets" % analysistype
    # Run the baiting process
    baitrprocesses(analysistype)
    print "\nIndexing %s targets" % analysistype
    # Index the combined target file
    SMALT.smaltindextargetsprocesses(armidatabase, currenttargetpath)
    print '\nPerforming %s reference mapping' % analysistype
    # Use SMALT to perform reference mapping of the combined target file
    SMALT.smaltmappingprocesses(seqdict, analysistype, 'SMALT')
    print "\nSorting mapped %s files" % analysistype
    # Use samtools to sort the bam file
    bamProcessor.sortingprocesses(seqdict, analysistype)
    print '\nIndexing sorted %s files' % analysistype
    # Use samtools to index the sorted bam file
    bamProcessor.bamindexingprocesses(seqdict, analysistype)
    print '\nParsing %s results' % analysistype
    # Use pysamstats to parse the bam files. Mike's armi module is called from within bamPysamStatsCombined
    parseddict = bamPysamStats.bamparseprocesses(seqdict, analysistype)
    bamPysamStatsCombined.armiparser(parseddict, seqdict, analysistype, reportfolder)
    # print json.dumps(seqdict, sort_keys=True, indent=4, separators=(',', ': '))


def rmlst():
    """Performs the necessary analyses on strains using armi targets"""
    global targetpath, seqdict
    # Set the analysis type
    analysistype = "rMLST"
    # Set the profile name
    profiletype = "rMLSTprofile"
    # Set the path of the analysistype data
    currenttargetpath = "%s%s" % (targetpath, analysistype)
    # Set the allele and profile path variables
    rmlstpath = currenttargetpath + "/alleles"
    rmlstprofilepath = currenttargetpath + "/profile"
    # Find the .json version of the profile file
    profilefile = glob("%s/*.json" % rmlstprofilepath)
    # If there is no .json file, use the .txt profile file instead. Note that this file must be edited (by removing
    # any columns after the last gene)
    if not profilefile:
        profilefile = glob("%s/*.txt" % rmlstprofilepath)[0]
    else:
        profilefile = profilefile[0]
    # Because the combined database of alleles is so large (~143 MB), it is filtered with uSearch prior to baiting
    # this is not done automatically, so the set-up of the bait files is slightly different here
    # In order to save time, a precomputed hash file is used
    baitfile = glob("%s/bait/*.gz" % rmlstpath)[0]
    baittargets(rmlstpath, analysistype)
    # Get a list of the targets into a variable - use a list comprehension to remove any database files present
    # target.fa, or target.fasta should be acceptable, but target.fa.nsi, or target.fasta.smi should be eliminated
    targets = [fasta for fasta in glob("%s/*.fa*" % rmlstpath) if not re.search("\.fa.*[\.]", fasta) and
               "Concatenated" not in fasta]
    # Add the bait file, the profile file, and a list of the targets to seqdict
    for strain in seqdict:
        seqdict[strain]["bait"]["fastqFiles"][analysistype] = baitfile
        seqdict[strain]["targets"][analysistype] = targets
        seqdict[strain][profiletype] = profilefile
    # Index the targets
    print "Indexing %s targets" % analysistype
    SMALT.smaltindextargetsprocesses(targets, rmlstpath)
    # Bait!
    print "\nFiltering .fastq files with %s targets" % analysistype
    baitrprocesses(analysistype)
    # Get the MLST profile into a dictionary
    print "\nLoading %s profiles" % analysistype
    profiledict, genedict = rawMLST.profilR(seqdict, profiletype)
    # Use SMALT to perform reference mapping
    print '\nPerforming %s reference mapping' % analysistype
    SMALT.smaltmappingprocesses(seqdict, analysistype, "SMALT")
    # Use samtools to sort bam files
    print "\nSorting mapped %s files" % analysistype
    bamProcessor.sortingprocesses(seqdict, analysistype)
    # Use samtools to index bam files
    print '\nIndexing sorted %s files' % analysistype
    bamProcessor.bamindexingprocesses(seqdict, analysistype)
    # Use pysamstats to parse bam files
    print '\nParsing %s results' % analysistype
    mlstmatches = bamPysamStats.bamparseprocesses(seqdict, analysistype)
    # Determine sequence types
    print '\nFinding multilocus sequence types'
    sequencetypes = rawMLST.sequenceTyper(mlstmatches, profiledict, genedict)
    # Create a report
    rawMLST.rMLSTreportMaker(seqdict, sequencetypes, analysistype, reportfolder, path)


def customtargets():
    """Performs genesipping on user-provided custom targets"""
    global targetpath, seqdict, args
    # Set the analysis type variable to 16S. This variable is important for retrieving 16S-specific data from seqdict
    analysistype = "custom"
    # Set the path of the analysistype data
    currenttargetpath = args['customTargetPath']
    # Create the bait target file (if necessary)
    baittargets(currenttargetpath, analysistype)
    # In order to save time, during baiting, a precomputed hash file is used
    hashfile = glob("%s/bait/*.gz" % currenttargetpath)
    # If the hashfile is present, use it
    if hashfile:
        baitfile = hashfile[0]
    # Otherwise, use the concatenated bait file
    else:
        baitfile = glob("%s/bait/*.fa*" % currenttargetpath)[0]
    print "Filtering .fastq files with %s targets" % analysistype
    # Get a list of the targets into a variable - use a list comprehension to remove any database files present
    # target.fa, or target.fasta should be acceptable, but target.fa.nsi, or target.fasta.smi should be eliminated
    customdatabase = [fasta for fasta in glob("%s/*.fa*" % currenttargetpath) if not re.search("\.fa.*[\.]", fasta) and
                      "Concatenated" not in fasta]
    # Populate seqdict with the name and path of the bait file used as well as the database
    for foldername in seqdict:
        seqdict[foldername]["bait"]["fastqFiles"][analysistype] = baitfile
        seqdict[foldername]["targets"][analysistype] = customdatabase
    # Run the baiting process
    baitrprocesses(analysistype)
    # Perform SMALT indexing of targets
    print "\nIndexing %s targets" % analysistype
    SMALT.smaltindextargetsprocesses(customdatabase, currenttargetpath)
    # Perform reference mapping with SMALT
    print '\nPerforming %s reference mapping' % analysistype
    SMALT.smaltmappingprocesses(seqdict, analysistype, "SMALT")
    # Use samtools to sort the bam files
    print "\nSorting mapped %s files" % analysistype
    bamProcessor.sortingprocesses(seqdict, analysistype)
    # Use samtools to index the sorted bam files
    print '\nIndexing sorted %s files' % analysistype
    bamProcessor.bamindexingprocesses(seqdict, analysistype)
    # Use pysamstats to parse the sorted, indexed bam files
    print '\nParsing %s results' % analysistype
    customdict = bamPysamStats.bamparseprocesses(seqdict, analysistype)
    # Create reports of the results
    reportr(customdict, analysistype)
    # import json
    # print json.dumps(customdict, sort_keys=True, indent=4, separators=(',', ': '))


def runner():
    """Calls the geneSipping functions in the appropriate order"""
    global sequencepath, path, miseqpath, projectname, forwardreads, reversereads, args
    # Initialise dictionary and list to store 16S typing results
    organismdict = defaultdict(make_dict)
    organismlist = []
    # Run the fastqCreator function if necessary
    if miseqpath:
        printtime("fastqCreator")
        fastqCreator.createfastq(miseqpath, miseqfolder, path, projectname,
                                 forwardreads, reversereads, customsamplesheet)
    # Run name extractor to determine the name of the samples. Additionally, this function decompresses any .gz files
    print "Finding sample names"
    nameextractor()
    # Run the 16S typing function if one or more of a few arguments are provided
    if args['16Styping'] or args['Mlst'] or args['pathotYping'] or args['Serotyping'] or args['Virulencetyping']:
        printtime("16S")
        organismdict, organismlist = sixteens()
    # Perform specified analyses
    if args['Mlst']:
        # Run MLST analyses
        printtime("MLST")
        mlst(organismdict, organismlist)
    if args['pathotYping']:
        # Perform pathotyping analyses
        printtime("pathotype")
        pathotyper(organismdict, organismlist, "pathotype")
    if args['Serotyping']:
        # Perform serotyping analyses
        printtime("serotype")
        pathotyper(organismdict, organismlist, "serotype")
    if args['Virulencetyping']:
        # Perform virulence typing
        printtime("virulence typing")
        virulencetyper(organismdict, "virulencetype")
    if args['armi']:
        # Perform armi analyses
        printtime("ARMI")
        armi()
    if args['rmlst']:
        # Perform rMLST analyses
        printtime("rMLST")
        rmlst()
    if args['customTargetPath']:
        # Analyse the custom targets
        printtime("customTargets")
        customtargets()

# Define the start time
start = time.time()

# Run the script
runner()

# Print a bold, green exit statement
print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'
# print json.dumps(seqdict, sort_keys=True, indent=4, separators=(',', ': '))

import pysamstats
import os
import sys
import json
import errno
import time
import subprocess
from collections import defaultdict
from glob import glob
from multiprocessing import Pool
from ARMICARD import decipher

__author__ = 'adamkoziol'


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    try:
        os.makedirs(inpath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


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


# def filler(listofdictionaries):
#     """
#     Properly populates the dictionary - when I tried to populate the dictionary within the multi-processed functions,
#     it was returned as a list (likely due to how the files were returned. This essentially
#      iterates through the list, and updates a dictionary appropriately
#      :param listofdictionaries: list of dictionaries returned from a multiprocessed function
#      """
#     # Initialise the dictionary
#     replacementdictionary = defaultdict(make_dict)
#     # Iterate through listofdictionaries
#     for dictionary in listofdictionaries:
#         # Update the dictionary with the dictionaries stored in the list
#         replacementdictionary.update(dictionary)
#     # Return the beautifully-populated dictionary
#     return replacementdictionary

def filler(listofdictionaries):
    """
    Properly populates the dictionary - when I tried to populate the dictionary within the multi-processed functions,
    it was returned as a list (likely due to how the files were returned. This essentially
     iterates through the list, and updates a dictionary appropriately
     :param listofdictionaries: list of dictionaries returned from a multiprocessed function
     """
    # Initialise the dictionary
    replacementdictionary = defaultdict(make_dict)
    # Iterate through listofdictionaries
    for dictionary in listofdictionaries:
        # Each dictionary follows a strict formula - iterate through all the nested values
        for strain in dictionary:
            for target in dictionary[strain]:
                for allele in dictionary[strain][target]:
                    for pos in dictionary[strain][target][allele]:
                        for depth, quality in dictionary[strain][target][allele][pos].iteritems():
                            # Populate the replacement dictionary with the values
                            replacementdictionary[strain][target][allele][pos][depth] = quality
    # Return the beautifully-populated dictionary
    return replacementdictionary


# Initialise parsedict
parsedict = defaultdict(make_dict)


def bestmatch(parseddict, seqdict):
    """
    Determines the best match to a target based on sequence identity
    :param parseddict: dictionary of filtered, parsed results
    :param seqdict: dictionary containing names and paths of files and folders required for these analyses
    """
    # Initialise a dictionary to hold the results
    resultdict = defaultdict(make_dict)
    # Iterate through the strains in the analysis
    for strain in seqdict:
        # Get the bait type value from seqdict
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        # Determine the location of the baited fastq files from seqdict
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        # Define the name of the file to hold the parsed results - allows for much quicker re-analysis of the data
        jsonprofile = "%s/%s_%s.json" % (fastqdir, strain, baittype)
        # If the file doesn't exist, perform the parsing
        if not os.path.isfile(jsonprofile):
            # Create a dictionary that will not store the strain name - this is useful in later parsing steps
            tempdict = defaultdict(make_dict)
            # Analyse each target
            for targetname in parseddict[strain]:
                # Initialise variables at 0. These variables will eventually store the best percent identity and depth
                # values for each target examined
                maxidentity = 0
                maxdepth = 0
                # Iterate through the records in parseddict
                for record in parseddict[strain][targetname]:
                    # Iterate through depth values
                    for depth, identity in parseddict[strain][targetname][record].iteritems():
                        # Iterate through identity values
                        # for identity in parseddict[strain][targetname][record][depth]:
                        # Find the highest identity value in the dictionary
                        if identity > maxidentity:
                            # If the current identity is greater than the highest observed identity,
                            # set the highest observed to the current identity
                            maxidentity = identity
                # Now that the highest identity value is known, iterate through the records again
                for record in parseddict[strain][targetname]:
                    # Iterate through the depth values
                    for depth, identity in parseddict[strain][targetname][record].iteritems():
                        # for depth in parseddict[strain][targetname][record]:
                        # Iterate through the identity values
                        # for identity in parseddict[strain][targetname][record][depth]:
                        # If the identity is equal to the largest identity
                        if identity == maxidentity:
                            # Find the best depth value
                            if depth > maxdepth:
                                # Set the maximum depth to the current depth
                                maxdepth = depth
                # Now that the highest identity and depth values are known, iterate through the records again
                for record in parseddict[strain][targetname]:
                    # Iterate through the depth values
                    for depth, identity in parseddict[strain][targetname][record].iteritems():
                        # for depth in parseddict[strain][targetname][record]:
                        # Iterate through the identity values
                        # for identity in parseddict[strain][targetname][record][depth]:
                        # If this records has the highest identity, and the highest depth value observed for that
                        # identity, then add the appropriate information to dictionaries
                        if depth == maxdepth:
                            # resultDict will store this best result
                            resultdict[strain][targetname][record][identity] = depth
                            # tempDict stores the same data as resultDict except for the strain name
                            tempdict[targetname][record][identity] = depth
            # Write the data to file
            jsonreport = open(jsonprofile, "wb")
            # Create the JSON-formatted output to write to file
            output = json.dumps(tempdict, sort_keys=True, indent=4, separators=(',', ': '))
            jsonreport.write(output)
            jsonreport.close()
        # If the JSON file already exists
        else:
            # Open the JSON file
            with open(jsonprofile, "rb") as jsonreport:
                # Load the data
                resultdict[strain].update(json.load(jsonreport))
    # Return the dictionary
    return resultdict


def armiparser(parseddict, seqdict, analysistype, reportfolder):
    """Creates a dictionary that can be put into Mike's ARMI decipher function
    :param parseddict: dictionary of parsed results
    :param seqdict: dictionary containing import paths and file names
    :param analysistype: string of the current analysis type
    :param reportfolder: folder in which reports are stored
    """
    # Initialise the dictionary
    targetdict = defaultdict(make_dict)
    # Iterate through the strains in the analysis
    for strain in seqdict:
        # Get the identity cutoff from seqdict - need to get baittype first
        identitycutoff = seqdict[strain]["cutoff"][analysistype]
        # Iterate through all the targets
        for target in seqdict[strain]["targets"][analysistype]:
            # Initialise targetpresent as false - this will change only if the percent identity is greater than the
            # identity cutoff
            targetpresent = False
            # Create the targetname variable from the target
            targetname = os.path.basename(target).split(".")[0]
            # Iterate through all the alleles for each target in parseddict
            for allele in parseddict[strain][target]:
                # Initialise the totaldepth and the number of nonsnps (number of matches to the reference)
                totaldepth = 0
                nonsnps = 0
                # Retrieve the length of the allele
                contiglength = len(seqdict[strain]["targetSequences"][analysistype][target]["allele"][allele])
                # Iterate through each position in the allele
                for pos in parseddict[strain][target][allele]:
                    # Number of matches
                    for matches in parseddict[strain][target][allele][pos]:
                        # Number of mismatches and depth of coverage
                        for mismatches, depth in parseddict[strain][target][allele][pos][matches].iteritems():
                            # Each position represents a non-SNP due to pre-filtering
                            nonsnps += 1
                            # Increment the total depth
                            totaldepth += depth
                # Calculate the total percent identity
                currentidentity = float("%.2f" % (float(nonsnps)/contiglength * 100))
                # If this identity is greater than the cutoff
                if currentidentity >= identitycutoff:
                    # The target is present in the strain
                    targetpresent = True
            # If the target is present, add a plus to Dict
            if targetpresent:
                targetdict[strain][targetname] = ["+"]
    # Set the path of the resistance dictionary
    antidict = json.load(open("/media/nas0/Jackson/ARMI_Docker/ARMI/aro3.json"))
    # Send the dictionaries, and report locations to the decipher function
    decipher(targetdict, antidict, reportfolder + "/geneSippr")
    return targetdict


def bamparseprocesses(seqdict, analysistype, reportfolder):
    """
    Multiprocessing for sorting bam files
    :param seqdict: dictionary containing import path and name information of files and folders
    :param analysistype: string of the analysis type
    :param reportfolder: string of output directory
    """
    # Initialise dictionary for later
    loadedresultsdict = defaultdict(make_dict)
    # Initialise argument list and pool for multiprocessing
    bamparseprocessesargs = []
    bamparseprocessespool = Pool()
    # Iterate through the strains for appending necessary arguments to list
    for strain in seqdict:
        # Set the identity cutoff depending on the analysis
        if analysistype == "pathotype":
            identitycutoff = 70
        elif analysistype == "MLST" or analysistype == "rMLST":
            identitycutoff = 95
        else:
            identitycutoff = 90
        # Retrieve necessary values from seqdict
        seqdict[strain]["cutoff"][analysistype] = identitycutoff
        # Retrieve bait type and determine the directory of the fastq files from seqdict
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        # Get the name of the target from the catfile variable
        catfile = seqdict[strain]["concatenatedTargets"][analysistype]
        targetname = os.path.basename(catfile).split(".")[0]
        # Set the input/output dir
        outdir = "%s/%s" % (fastqdir, targetname)
        # Make a list of the sorted bam files
        bamfiles = glob("%s/*_sorted.bam" % outdir)[0]
        # Set the name of the JSON file to store the results
        jsonprofile = "%s/%s_matchDict.json" % (fastqdir, analysistype)
        # If the analyses have not previously been performed, add the arguments to the appropriate list
        if not os.path.isfile(jsonprofile):
            bamparseprocessesargs.append((strain, catfile, bamfiles))
        # Otherwise, load the results from the json file
        else:
            with open(jsonprofile, "rb") as jsonreport:
                # Load the data
                loadedresultsdict[strain].update(json.load(jsonreport))
                dotter()
    # Run any necessary parsing
    parselist = bamparseprocessespool.map(bamparse, bamparseprocessesargs)
    # Create a dictionary from the returned list of dictionaries
    parseddict = filler(parselist)
    # Initialise dictionaries to be updated later
    resultsdict = defaultdict(make_dict)
    bestdict = defaultdict(make_dict)
    # Load the length of the targets using the .fai files generated in the bamParse function
    seqdict = targetlength(seqdict, analysistype)
    # Iterate through the strains in order to write the results to a JSON file
    for strain in seqdict:
        # Get the bait type
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        # Define the JSON profile file
        jsonprofile = "%s/%s_concatenatedMatchDict.json" % (fastqdir, analysistype)
        # If the file doesn't exist, create it, and fill it with results
        if not os.path.isfile(jsonprofile):
            jsonreport = open(jsonprofile, "wb")
            output = json.dumps(parseddict[strain], sort_keys=True, indent=4, separators=(',', ': '))
            jsonreport.write(output)
            jsonreport.close()
    # Update parseddict with any results loaded from file
    parseddict.update(loadedresultsdict)
    # ARMI analyses require special parsing compared to the other types of analysis
    if analysistype == "ARMI":
        armidict = armiparser(parseddict, seqdict, analysistype, reportfolder)
        return armidict
    else:
        # Parse and filter the results
        resultsdict.update(dictparser(parseddict, seqdict, analysistype))
        # Find the best result
        bestdict.update(bestmatch(resultsdict, seqdict))
        return bestdict


def bamparse((strain, target, bamfile)):
    """Parses bam files using pysam stats"""
    global parsedict
    # Use the stat_baseq_ext (extended base quality statistics) function of pysam stats to return records parsed
    # from sorted bam files
    for rec in pysamstats.stat_baseq_ext(alignmentfile=bamfile, fafile=target):
        # Values of interest can be retrieved using the appropriate keys
        # Simple filtering statement: if the number of matches at a particular position in the reference sequence is
        # greater than the number of mismatches, and the total depth is 5 or more, add the position of the results
        if rec['matches'] > rec['mismatches'] and rec['reads_all'] > 4:
            # Populate the dictionary with the appropriate values
            parsedict[strain][target][rec['chrom']][float(rec['pos'])][rec['reads_all']] = rec['rms_baseq']
    dotter()
    return parsedict


def dictparser(parseddict, seqdict, analysistype):
    """
    Filter the results that are below the identity cutoff threshold
    :param parseddict: dictionary containing filtered results
    :param seqdict: dictionary containing names and paths of important files and folders in the analysis
    :param analysistype: string of the analysis type (e.g. 16S)
    """
    # Imaginatively named dictionary to store the filtered results
    plusdict = defaultdict(make_dict)
    # Iterate through the strains with results
    for strain in parseddict:
        # Retrieve the identity cutoff value from seqDict
        identitycutoff = seqdict[strain]["cutoff"][analysistype]
        # Iterate through the targets
        for target in seqdict[strain]["targets"][analysistype]:
            # Set the target name
            targetname = os.path.basename(target).split(".")[0]
            # For every allele in the targets (one or more)
            # cattarget = parseddict[strain].keys()[0]
            # for allele in parseddict[strain][cattarget]:
            for allele in seqdict[strain]["targetSequences"][analysistype][target]["allele"]:
                # Initialise the variables to store the parsed values
                totaldepth = 0
                nonsnps = 0
                # Retrieve the contig (allele) length from seqDict
                contiglength = seqdict[strain]["targetSequences"][analysistype][target]["allele"][allele]
                # Iterate through the cat file
                # for cattarget in parseddict[strain]:
                cattarget = parseddict[strain].keys()[0]
                # Iterate through each individual position in the reference mapped bam files
                for pos in sorted(parseddict[strain][cattarget][allele]):
                    # Get the depth and quality values from the dictionary
                    for depth, quality in parseddict[strain][cattarget][allele][pos].iteritems():
                        # As these results were filtered in bamParse, every position represents a non-SNP
                        nonsnps += 1
                        # Increment the total depth value with the current depth
                        totaldepth += float(depth)
                # Calculate the average depth and identity values using the contig length
                averagedepth = float("%.2f" % (float(totaldepth)/contiglength))
                percentidentity = float("%.2f" % (float(nonsnps)/contiglength * 100))
                # If the observed percent identity is greater than the cutoff value, add the results to the dictionary
                if percentidentity >= identitycutoff:
                    plusdict[strain][targetname][allele][averagedepth] = percentidentity
    # Return the dictionary
    return plusdict


def targetlength(seqdict, analysistype):
    """
    Uses .fai files to determine the length of target sequences
    :param seqdict: dictionary containing names and paths of files and folders
    :param analysistype: string of the analysis type (e.g. 16S)
    """
    # Iterate through all the strains
    for strain in seqdict:
        # Iterate through the targets for each strain in seqDict
        for target in seqdict[strain]["targets"][analysistype]:
            # The fai file is the target will be the target name, with .fai appended on the end
            faifile = target + ".fai"
            # If the faidx file does not exist, create it
            if not os.path.isfile(faifile):
                faidxcommand = "samtools faidx %s" % target
                # Define /dev/null
                fnull = open(os.devnull, 'wb')
                # Run the command
                subprocess.call(faidxcommand, shell=True, stdout=fnull, stderr=fnull)
            # Open the .fai file
            with open(faifile, "rb") as faiInfo:
                # The .fai file is essentially a tab-delimited list of alleles and the length of the allele
                for line in faiInfo:
                    # Each line contains the data for one allele sequence - split on the tab
                    data = line.split("\t")
                    # Get the allele name (data[0]) and the length of the allele (data[1]) into seqDict
                    seqdict[strain]["targetSequences"][analysistype][target]["allele"][data[0]] = float(data[1])
        # Run dotter to show progress
        dotter()
    # Return the updated seqDict
    return seqdict

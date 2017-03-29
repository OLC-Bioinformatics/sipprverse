import pysamstats
import os
import sys
import json
import errno
import time
from collections import defaultdict
from glob import glob
from multiprocessing import Pool

__author__ = 'adamkoziol'


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)


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


def bestmatch(parseddict, seqdict, analysistype):
    """
    Determines the best match to a target based on sequence identity
    :param parseddict: dictionary of filtered, parsed results
    :param seqdict: dictionary containing names and paths of files and folders required for these analyses
    :param analysistype: string of the current analysis type
    """
    # Initialise a dictionary to hold the results
    resultdict = defaultdict(make_dict)
    # Iterate through the strains in the analysis
    for strain in seqdict:
        # Analyse each target
        for targetname in parseddict[strain]:
            # Initialise variables at 0. These variables will eventually store the highest percent identity and depth
            # values for each target examined
            maxidentity = 0
            maxdepth = 0
            # Get the bait type value from seqDict
            baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
            # Determine the location of the baited fastq files from seqDict
            fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
            # Define the name of the file to hold the parsed results - allows for much quicker re-analysis of the data
            jsonprofile = "%s/%s/%s_%s_%s.json" % (fastqdir, targetname, strain, baittype, targetname)
            # If the file doesn't exist, perform the parsing
            if not os.path.isfile(jsonprofile):
                # Create a dictionary that will not store the strain name - this is useful in later parsing steps
                tempdict = defaultdict(make_dict)
                # Iterate through the records in parsedDict
                for record in parseddict[strain][targetname]:
                    # Iterate through depth values
                    for depth in parseddict[strain][targetname][record]:
                        # Iterate through identity values
                        # Find the highest identity value in the dictionary
                        identity = parseddict[strain][targetname][record][depth]
                        if identity > maxidentity:
                            # If the current identity is greater than the highest observed identity,
                            # set the highest observed to the current identity
                            maxidentity = identity
                # Now that the highest identity value is known, iterate through the records again
                for record in parseddict[strain][targetname]:
                    # Iterate through the depth values
                    for depth in parseddict[strain][targetname][record]:
                        # Iterate through the identity values
                        # If the identity is equal to the largest identity
                        identity = parseddict[strain][targetname][record][depth]
                        # Virulence type is special - I want to allow multiple best matches and the matches don't
                        # have to all be the same - allowing a decrease of 5% from the best hit
                        if analysistype == "virulencetype":
                            if identity >= maxidentity - 5:
                                # resultDict will store this best result
                                resultdict[strain][targetname][record][identity] = depth
                                # tempDict stores the same data as resultDict except for the strain name
                                tempdict[targetname][record][identity] = depth
                        # Find the best depth value
                        else:
                            if identity == maxidentity:
                                # If the depth is equal to the best value of depth
                                if depth > maxdepth:
                                    # Set the maximum depth to the current depth
                                    maxdepth = depth
                # Now that the highest identity and depth values are known, iterate through the records again
                for record in parseddict[strain][targetname]:
                    # Iterate through the depth values
                    for depth in parseddict[strain][targetname][record]:
                        # Iterate through the identity values
                        identity = parseddict[strain][targetname][record][depth]
                        # If this records has the highest identity, and the highest depth value observed for that
                        # identity, then add the appropriate information to dictionaries
                        if depth == maxdepth and analysistype != "virulencetype":
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
                with open(jsonprofile, "rb") as jsonReport:
                    # Load the data
                    resultdict[strain].update(json.load(jsonReport))
    # Return the dictionary
    return resultdict


def classifyr(bestresult, seqdict, analysistype):
    """
    Determines genus of strains based on reference mapping results
    :param bestresult: dictionary containing the best match to a target
    :param seqdict: dictionary containing name and path information of important files and folders
    :param analysistype: string of analysis type
    """
    # Initialise variables
    genusdict = defaultdict(make_dict)
    genusset = set()
    # Iterate through strains
    for strain in seqdict:
        # Pull the path of the target files from seqDict
        targetpath = os.path.split(seqdict[strain]["targets"][analysistype][0])[0]
        # The taxa file must have a .tax file extension and be present in the target path
        taxafile = glob("%s/*.tax" % targetpath)[0]
        # Iterate through the targets in bestResult
        for targetname in bestresult[strain]:
            # Open the taxa file
            with open(taxafile) as taxa:
                for line in taxa:
                    # Looks to see if the accession number (with trailing 1 stripped off, and ".1" re-added to the end)
                    # of the best result is the same as the first entry of the line
                    if bestresult[strain][targetname].keys()[0][0:-1] + ".1" == line.split("\t")[0]:
                        # Set the genus as the second last value in the taxonomy string separated by ";"
                        # Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia;
                        # would give Escherichia. I don't remember why this value is also split on the "_"
                        genus = line.split(";")[-2].split("_")[0]
                        if genus == "Shigella":
                            genus = "Escherichia"
                        # Add the best results to the bestResult dictionary
                        for bestcontig in bestresult[strain]:
                            genusdict[strain][genus][bestcontig][bestresult[strain][bestcontig].keys()[0]] = \
                                bestresult[strain][bestcontig].values()[0]
                            genusset.add(genus)
    return genusdict, genusset


def bamparseprocesses(seqdict, analysistype):
    """
    Multiprocessing for parsing bam files
    :param seqdict: dictionary containing import path and name information of files and folders
    :param analysistype: string of the analysis type
    """
    # Initialise dictionary, argument list, and Pool
    loadedresultsdict = defaultdict(make_dict)
    bamparseprocessesargs = []
    bamparseprocessespool = Pool()
    # Iterate through the strains
    for strain in seqdict:
        # Define the cutoff values depending on the analysis. Default is 90 identity
        # Pathotyping analysis has a lower cutoff value, while (r)MLST has a higher value
        if analysistype == "pathotype":
            identitycutoff = 70
        elif analysistype == "MLST" or analysistype == "rMLST":
            identitycutoff = 95
        elif analysistype == 'virulencetype':
            identitycutoff = 98
        else:
            identitycutoff = 90
        #  Store the identityCutoff in seqDict
        seqdict[strain]["cutoff"][analysistype] = identitycutoff
        # Retrieve bait type and determine the directory of the fastq files from seqdict
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        # Set the name of the JSON file to store the results
        jsonprofile = "%s/%s_matchDict.json" % (fastqdir, analysistype)
        # If the JSON file hasn't been created, parse the bam files
        if not os.path.isfile(jsonprofile):
            for target in seqdict[strain]["targets"][analysistype]:
                # Get the name of the target from the target variable
                targetname = os.path.basename(target).split(".")[0]
                # Set the input/output dir
                outdir = "%s/%s" % (fastqdir, targetname)
                # Make a list of the sorted bam files
                try:
                    bamfile = glob("%s/*_sorted.bam" % outdir)[0]
                    # Append a tuple of the required arguments to the argument list
                    bamparseprocessesargs.append((strain, target, bamfile))
                except IndexError:
                    print outdir

        # If the JSON file exists, read the results from it rather than performing the parsing again
        else:
            # Open the JSON file
            with open(jsonprofile, "rb") as jsonreport:
                # Load the results from the JSON file into a dictionary
                loadedresultsdict[strain].update(json.load(jsonreport))
                dotter()
    # Run the multiprocessed bam parsing
    parselist = bamparseprocessespool.map(bamparse, bamparseprocessesargs)
    # Change the returned list of dictionaries into a nested dictionary
    parseddict = filler(parselist)
    # Load the length of the targets using the .fai files generated in the bamParse function
    seqdict = targetlength(seqdict, analysistype)
    # Iterate through the strains in order to write the results to a JSON file
    for strain in seqdict:
        # Get the bait type
        baittype = seqdict[strain]["bait"]["fastqFiles"].keys()[0]
        fastqdir = os.path.split(seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
        # Define the JSON profile file
        jsonprofile = "%s/%s_matchDict.json" % (fastqdir, analysistype)
        # If the file doesn't exist, create it, and fill it with results
        if not os.path.isfile(jsonprofile):
            jsonreport = open(jsonprofile, "wb")
            output = json.dumps(parseddict[strain], sort_keys=True, indent=4, separators=(',', ': '))
            jsonreport.write(output)
            jsonreport.close()
    # Update parseddict with any results loaded from file
    parseddict.update(loadedresultsdict)
    # Analyses are filtered based on whether results are above the identity cutoff
    filtereddict = dictparser(parseddict, seqdict, analysistype)
    # Find the best result
    resultsdict = bestmatch(filtereddict, seqdict, analysistype)
    # Find the genus from the best result for 16S/18S analyses
    if analysistype == "16S" or analysistype == "18S":
        genusdict, genusset = classifyr(resultsdict, seqdict, analysistype)
        # Return the dictionary with the genera, and a list of all genera encountered in these analyses
        return genusdict, list(genusset)
    # print json.dumps(filtereddict, sort_keys=True, indent=4, separators=(',', ': '))
    # Otherwise return resultsDict
    else:
        return resultsdict


def bamparse((strain, target, bamfile)):
    """Parses bam files using pysam stats"""
    parsedict = defaultdict(make_dict)
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
    # Imaginatively name dictionary to store the filtered results
    plusdict = defaultdict(make_dict)
    # Iterate through the strains with results
    for strain in parseddict:
        # Retrieve the identity cutoff value from seqDict
        identitycutoff = seqdict[strain]["cutoff"][analysistype]
        # Iterate through the targets
        for target in parseddict[strain]:
            # Set the target name
            targetname = os.path.basename(target).split(".")[0]
            # For every allele in the targets (one or more)
            for allele in parseddict[strain][target]:
                # Initialise the variables to store the parsed values
                totaldepth = 0
                nonsnps = 0
                # Retrieve the contig (allele) length from seqDict
                contiglength = seqdict[strain]["targetSequences"][analysistype][target]["allele"][allele]
                # Iterate through each individual position in the reference mapped bam files
                for pos in sorted(parseddict[strain][target][allele]):
                    # Get the depth and quality values from the dictionary
                    for depth, quality in parseddict[strain][target][allele][pos].iteritems():
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
    import subprocess
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

#!/usr/bin/env python
from collections import defaultdict
from multiprocessing import Pool
import bamPysamStats
import os
import json


__author__ = 'adamkoziol'


class ArmiParse(object):

    def parsearmi(self):
        """
        Multiprocessing for parsing bam files
        :param seqdict: dictionary containing import path and name information of files and folders
        :param analysistype: string of the analysis type
        """
        from glob import glob
        # Initialise dictionary, argument list, and Pool
        loadedresultsdict = defaultdict(bamPysamStats.make_dict)
        bamparseprocessesargs = []
        bamparseprocessespool = Pool()
        # Iterate through the strains
        for strain in self.seqdict:

            #  Store the identityCutoff in seqDict
            self.seqdict[strain]["cutoff"][self.analysistype] = self.identitycutoff
            # Retrieve bait type and determine the directory of the fastq files from seqdict
            baittype = self.seqdict[strain]["bait"]["fastqFiles"].keys()[0]
            fastqdir = os.path.split(self.seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
            # Set the name of the JSON file to store the results
            jsonprofile = "%s/%s_matchDict.json" % (fastqdir, self.analysistype)
            # If the JSON file hasn't been created, parse the bam files
            if not os.path.isfile(jsonprofile):
                bamfile = glob('{}/ARMIConcatenated/*_sorted.bam'.format(fastqdir))[0]
                for target in self.seqdict[strain]["targets"][self.analysistype]:
                    # Get the name of the target from the target variable
                    targetname = os.path.basename(target).split(".")[0]
                    # Set the input/output dir
                    # outdir = "%s/%s" % (fastqdir, targetname)
                    # Make a list of the sorted bam files
                    # bamfile = glob("%s/*_sorted.bam" % outdir)[0]
                    # catfile = self.seqdict[strain]["concatenatedTargets"][self.analysistype]
                    # print "!", catfile


                    # Append a tuple of the required arguments to the argument list
                    bamparseprocessesargs.append((strain, target, bamfile, targetname))
                    # print strain, target, bamfile
            # If the JSON file exists, read the results from it rather than performing the parsing again
            else:
                # Open the JSON file
                with open(jsonprofile, "rb") as jsonreport:
                    # Load the results from the JSON file into a dictionary
                    loadedresultsdict[strain].update(json.load(jsonreport))
                    bamPysamStats.dotter()
        # Run the multiprocessed bam parsing
        parselist = bamparseprocessespool.map(bamparse, bamparseprocessesargs)
        # Change the returned list of dictionaries into a nested dictionary
        self.parseddict = bamPysamStats.filler(parselist)
        # Load the length of the targets using the .fai files generated in the bamParse function
        self.seqdict = bamPysamStats.targetlength(self.seqdict, self.analysistype)
        # Iterate through the strains in order to write the results to a JSON file
        for strain in self.seqdict:
            # Get the bait type
            baittype = self.seqdict[strain]["bait"]["fastqFiles"].keys()[0]
            fastqdir = os.path.split(self.seqdict[strain]["bait"]["fastqFiles"][baittype][0])[0]
            # Define the JSON profile file
            jsonprofile = "%s/%s_matchDict.json" % (fastqdir, self.analysistype)
            # If the file doesn't exist, create it, and fill it with results
            if not os.path.isfile(jsonprofile):
                jsonreport = open(jsonprofile, "wb")
                output = json.dumps(self.parseddict[strain], sort_keys=True, indent=4, separators=(',', ': '))
                jsonreport.write(output)
                jsonreport.close()
        self.armiparser()
        # print json.dumps(parselist, sort_keys=True, indent=4, separators=(',', ': '))

    def armiparser(self):
        """Creates a dictionary that can be put into Mike's ARMI decipher function
        :param parseddict: dictionary of parsed results
        :param seqdict: dictionary containing import paths and file names
        :param analysistype: string of the current analysis type
        :param reportfolder: folder in which reports are stored
        """
        from ARMICARD import decipher
        # Initialise a variable to store the target path used in creating the dictionary of antimicrobial resistances
        targetpath = ""
        # Initialise the dictionary
        targetdict = defaultdict(bamPysamStats.make_dict)
        # Iterate through the strains in the analysis
        for strain in self.seqdict:
            # Iterate through all the targets
            for target in sorted(self.seqdict[strain]["targets"][self.analysistype]):
                # Initialise targetpresent as false - this will change only if the percent identity is greater than the
                # identity cutoff
                targetpresent = False
                # Create the targetname variable from the target
                targetname = os.path.basename(target).split(".")[0]
                targetpath = os.path.split(target)[0]
                # Iterate through all the alleles for each target in parseddict
                for allele in self.parseddict[strain][target]:
                    # Initialise the totaldepth and the number of nonsnps (number of matches to the reference)
                    totaldepth = 0
                    nonsnps = 0
                    # Retrieve the length of the allele
                    contiglength = len(self.seqdict[strain]["targetSequences"][self.analysistype]
                                       [target]["allele"][allele])
                    # Iterate through each position in the allele
                    for pos in self.parseddict[strain][target][allele]:
                        # Number of matches
                        for matches in self.parseddict[strain][target][allele][pos]:
                            # Number of mismatches and depth of coverage
                            for mismatches, depth in self.parseddict[strain][target][allele][pos][matches].iteritems():
                                # Each position represents a non-SNP due to pre-filtering
                                nonsnps += 1
                                # Increment the total depth
                                totaldepth += depth
                    # Calculate the total percent identity
                    currentidentity = float("%.2f" % (float(nonsnps)/contiglength * 100))
                    print allele, currentidentity
                    # If this identity is greater than the cutoff
                    if currentidentity >= self.identitycutoff:
                        # The target is present in the strain
                        targetpresent = True
                # If the target is present, add a plus to Dict
                if targetpresent:
                    targetdict[strain][targetname] = ["+"]
        # Set the path of the resistance dictionary
        # antidict = json.load(open("/media/nas0/Jackson/ARMI_Docker/ARMI/aro3.json"))
        antidict = json.load(open("%s/aro3.json" % targetpath))
        # Send the dictionaries, and report locations to the decipher function
        # import json
        # print json.dumps(targetdict, sort_keys=True, indent=4, separators=(',', ': '))
        decipher(targetdict, antidict, self.reportfolder + "/geneSippr")
        return targetdict

    def __init__(self, seqdict, analysistype, reportfolder):
        self.seqdict = seqdict
        self.analysistype = analysistype
        self.reportfolder = reportfolder
        self.identitycutoff = 90
        self.parseddict = bamPysamStats.defaultdict(bamPysamStats.make_dict)
        self.parsearmi()


def bamparse((strain, target, bamfile, targetname)):
    """Parses bam files using pysam stats"""
    import pysamstats
    parsedict = defaultdict(bamPysamStats.make_dict)
    # Use the stat_baseq_ext (extended base quality statistics) function of pysam stats to return records parsed
    # from sorted bam files
    for rec in pysamstats.stat_baseq_ext(alignmentfile=bamfile, fafile=target):
        # Values of interest can be retrieved using the appropriate keys
        # Simple filtering statement: if the number of matches at a particular position in the reference sequence is
        # greater than the number of mismatches, and the total depth is 5 or more, add the position of the results
        if rec['matches'] > rec['mismatches'] and rec['reads_all'] > 4:
            # Populate the dictionary with the appropriate values
            parsedict[strain][target][rec['chrom']][float(rec['pos'])][rec['reads_all']] = rec['rms_baseq']
    bamPysamStats.dotter()

    # dotter()
    return parsedict

#!/usr/bin/env python
# Import the necessary modules
from SPAdesPipeline.OLCspades.accessoryFunctions import *
import SPAdesPipeline.OLCspades.metadataprinter as metadataprinter
import allelefind
import customtargets
import versions

__author__ = 'adamkoziol'


class Vtyper(object):
    def vtyper(self):
        """
        Calls the necessary methods in the appropriate order
        """
        import vtyperesults
        import createObject
        # Create a sample object
        self.runmetadata = createObject.ObjectCreation(self)
        # Run the baiting, mapping, sorting, and parsing method. Include a cutoff of 1.0 and a match bonus of 5
        custom = customtargets.Custom(self, 'vtyper', 1.0, 5)
        custom.targets()
        # Create the report
        vtyperesults.VtypeResults(self, 'vtyper')
        # Get the versions of the software used
        versions.Versions(self)
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: command line arguments
        :param pipelinecommit: git commit
        :param startingtime: starting time of the analyses
        :param scriptpath: the home path of the script
        """
        import multiprocessing
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.starttime = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        self.args = args
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.sequencepath = os.path.join(args.sequencepath, '')
        assert os.path.isdir(self.sequencepath), u'Supplied sequence path is not a valid directory {0!r:s}' \
            .format(self.sequencepath)
        self.targetpath = os.path.join(args.targetpath, '')
        assert os.path.isdir(self.targetpath), u'Supplied target path is not a valid directory {0!r:s}' \
            .format(self.targetpath)
        # In order to reuse code, the custom target path, which is used in the customtargets method must be used
        self.customtargetpath = self.targetpath
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = args.threads if args.threads else multiprocessing.cpu_count()
        # Determine if the determining the alleles is necessary
        self.alleles = args.alleles
        # Initialise the metadata
        self.runmetadata = MetadataObject()
        # Run the analyses
        self.vtyper()
        # Optionally run the allele-finding methods
        if self.alleles:
            alleles = allelefind.FindAllele(self, 'vtalleles')
            allelefind.FindAllele.allelefinder(alleles)
            # Print the metadata to file
            metadataprinter.MetadataPrinter(self)

if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    import subprocess
    import time

    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Parser for arguments
    parser = ArgumentParser(description='Verotoxin subtyping on raw reads')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    parser.add_argument('-t', '--targetpath',
                        required=True,
                        help='Path of target files to process.')
    parser.add_argument('-n', '--threads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-d', '--detailedReports',
                        action='store_true',
                        help='Provide detailed reports with percent identity and depth of coverage values '
                             'rather than just "+" for positive results')
    parser.add_argument('-a', '--alleles',
                        action='store_true',
                        help='Determine which vtx allele(s) are present in each strain. Save new alleles to be curated')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Define the start time
    start = time.time()
    # Run the script
    Vtyper(arguments, commit, start, homepath)
    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'

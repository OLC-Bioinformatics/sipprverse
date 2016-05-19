#!/usr/bin/env python
# Import the necessary modules
from SPAdesPipeline.OLCspades.accessoryFunctions import *

__author__ = 'adamkoziol'


class Vtyper(object):
    def vtyper(self):
        """
        Calls the necessary methods in the appropriate order
        """
        import customtargets
        import vtyperesults
        import createObject
        import SPAdesPipeline.OLCspades.metadataprinter as metadataprinter
        # Create a sample object
        self.runmetadata = createObject.ObjectCreation(self)
        # Run the baiting, mapping, sorting, and parsing method. Include a cutoff of 1.0 and a match bonus of 5
        customtargets.Custom(self, 'vtyper', 1.0, 5)
        # Create the report
        vtyperesults.VtypeResults(self, 'vtyper')
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
        self.runmetadata = MetadataObject()
        # Run the analyses
        self.vtyper()


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

    # Get the arguments into an object
    arguments = parser.parse_args()
    # Define the start time
    start = time.time()
    # Run the script
    Vtyper(arguments, commit, start, homepath)
    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'

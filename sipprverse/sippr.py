#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, printtime
from accessoryFunctions.metadataprinter import MetadataPrinter
from sixteenS.sixteens_full import SixteenS as SixteensFull
from sipprCommon.objectprep import Objectprep
from sipprCommon.sippingmethods import Sippr
from serosippr.serosippr import SeroSippr
from reporter.reports import Reports
from argparse import ArgumentParser
import multiprocessing
import subprocess
import time
import os

__author__ = 'adamkoziol'


class Sipprverse(object):

    def main(self):
        """
        Run the necessary methods in the correct order
        """
        printtime('Starting {} analysis pipeline'.format(self.analysistype), self.starttime)
        # Create the objects to be used in the analyses
        objects = Objectprep(self)
        objects.objectprep()
        self.runmetadata = objects.samples
        self.threads = int(self.cpus / len(self.runmetadata.samples)) if self.cpus / len(self.runmetadata.samples) > 1 \
            else 1
        # Run the genesippr analyses
        self.analysistype = 'genesippr'
        self.targetpath = os.path.join(self.reffilepath, self.analysistype, '')
        Sippr(self, self.cutoff)
        # Create the reports
        self.reports = Reports(self)
        Reports.reporter(self.reports)
        # Run the 16S analyses using the filtered database
        self.targetpath = self.reffilepath
        # Run the 16S analyses
        self.analysistype = 'sixteens_full'
        SixteensFull(self, self.commit, self.starttime, self.homepath, 'sixteens_full', 0.985)
        # Run the GDCS analysis
        self.analysistype = 'GDCS'
        self.pipeline = True
        Sippr(self, 0.95)
        # Create the reports
        Reports.gdcsreporter(self.reports)
        # Perform serotyping for samples classified as Escherichia
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                sample.mash = GenObject()
                try:
                    sample.mash.closestrefseqgenus = sample.general.closestrefseqgenus
                    for genus, species in self.taxonomy.items():
                        if genus == sample.mash.closestrefseqgenus:
                            sample.mash.closestrefseqspecies = species
                except KeyError:
                    sample.mash.closestrefseqgenus = 'NA'
                    sample.mash.closestrefseqspecies = 'NA'
            else:
                sample.mash.closestrefseqgenus = 'NA'
                sample.mash.closestrefseqspecies = 'NA'
        SeroSippr(self, self.commit, self.starttime, self.homepath, 'serosippr', 0.95, True)
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: command line arguments
        :param pipelinecommit: pipeline commit or version
        :param startingtime: time the script was started
        :param scriptpath: home path of the script
        """
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.starttime = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.sequencepath = os.path.join(args.sequencepath, '')
        self.seqpath = self.sequencepath
        self.targetpath = os.path.join(args.targetpath, '')
        # ref file path is used to work with sub module code with a different naming scheme
        self.reffilepath = self.targetpath
        self.reportpath = os.path.join(self.path, 'reports')
        make_path(self.reportpath)
        assert os.path.isdir(self.targetpath), u'Target path is not a valid directory {0!r:s}' \
            .format(self.targetpath)
        self.bcltofastq = args.bcl2fastq
        self.miseqpath = args.miseqpath
        self.miseqfolder = args.miseqfolder
        self.fastqdestination = args.destinationfastq
        self.forwardlength = args.readlengthforward
        self.reverselength = args.readlengthreverse
        self.numreads = 2 if self.reverselength != 0 else 1
        self.customsamplesheet = args.customsamplesheet
        # Set the custom cutoff value
        self.cutoff = float()
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.numthreads if args.numthreads else multiprocessing.cpu_count())
        self.threads = int()
        self.runmetadata = MetadataObject()
        self.taxonomy = {'Escherichia': 'coli', 'Listeria': 'monocytogenes', 'Salmonella': 'enterica'}
        self.analysistype = 'GeneSippr'
        self.copy = args.copy
        self.pipeline = False
        self.forward = str()
        self.reverse = str()
        self.index = str()
        self.header = dict()
        self.rundata = dict()
        self.completed = list()
        self.incomplete = list()
        self.analysescomplete = False
        self.final = False
        self.sum = int()
        self.completemetadata = list()
        self.samplesheetpath = str()
        self.samples = list()
        self.logfile = os.path.join(self.path, 'log')
        self.reports = str()
        # Run the method
        self.main()


if __name__ == '__main__':
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Parser for arguments
    parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    parser.add_argument('-t', '--targetpath',
                        required=True,
                        help='Path of target files to process.')
    parser.add_argument('-n', '--numthreads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-b', '--bcl2fastq',
                        action='store_true',
                        help='Optionally run bcl2fastq on an in-progress Illumina MiSeq run. Must include:'
                             'miseqpath, and miseqfolder arguments, and optionally readlengthforward, '
                             'readlengthreverse, and projectName arguments.')
    parser.add_argument('-m', '--miseqpath',
                        help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', '--miseqfolder',
                        help='Name of the folder containing MiSeq run data')
    parser.add_argument('-d', '--destinationfastq',
                        help='Optional folder path to store .fastq files created using the fastqCreation module. '
                             'Defaults to path/miseqfolder')
    parser.add_argument('-r1', '--readlengthforward',
                        default='full',
                        help='Length of forward reads to use. Can specify "full" to take the full length of '
                             'forward reads specified on the SampleSheet')
    parser.add_argument('-r2', '--readlengthreverse',
                        default='full',
                        help='Length of reverse reads to use. Can specify "full" to take the full length of '
                             'reverse reads specified on the SampleSheet')
    parser.add_argument('-c', '--customsamplesheet',
                        help='Path of folder containing a custom sample sheet (still must be named "SampleSheet.csv")')
    parser.add_argument('-P', '--projectName',
                        help='A name for the analyses. If nothing is provided, then the "Sample_Project" field '
                             'in the provided sample sheet will be used. Please note that bcl2fastq creates '
                             'subfolders using the project name, so if multiple names are provided, the results '
                             'will be split as into multiple projects')
    parser.add_argument('-D', '--detailedReports',
                        action='store_true',
                        help='Provide detailed reports with percent identity and depth of coverage values '
                             'rather than just "+" for positive results')
    parser.add_argument('-u', '--customcutoffs',
                        default=0.8,
                        help='Custom cutoff values')
    parser.add_argument('-C', '--copy',
                        action='store_true',
                        help='Normally, the program will create symbolic links of the files into the sequence path, '
                             'however, the are occasions when it is necessary to copy the files instead')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    Sipprverse(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')

#!/usr/bin/env python
import subprocess
from sipprCommon.sippingmethods import Sippr
from sipprCommon.objectprep import Objectprep
from accessoryFunctions.accessoryFunctions import printtime, make_path, MetadataObject
from accessoryFunctions.metadataprinter import MetadataPrinter
# Argument parser for user-inputted values, and a nifty help menu
from argparse import ArgumentParser
import multiprocessing
import os
import time
__author__ = 'adamkoziol'


class GeneSippr(object):

    def runner(self):
        """
        Run the necessary methods in the correct order
        """
        printtime('Starting {} analysis pipeline'.format(self.analysistype), self.starttime)
        if not self.pipeline:
            general = None
            for sample in self.runmetadata.samples:
                general = getattr(sample, 'general')
            if general is None:
                # Create the objects to be used in the analyses
                objects = Objectprep(self)
                objects.objectprep()
                self.runmetadata = objects.samples
        # Run the analyses
        Sippr(self, self.cutoff)
        # Create the reports
        self.reporter()
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create the path in which the reports are stored
        printtime('Creating reports', self.starttime)
        make_path(self.reportpath)
        data = 'Strain,Gene,PercentIdentity,FoldCoverage\n'
        with open(os.path.join(self.reportpath, '{}.csv'.format(self.analysistype)), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                try:
                    if sample[self.analysistype].results:
                        multiple = False
                        for name, identity in sample[self.analysistype].results.items():
                            if multiple:
                                data += ','
                            data += '{},{},{}\n'.format(name, identity, sample[self.analysistype].avgdepth[name])
                            multiple = True
                    else:
                        data += '\n'
                except KeyError:
                    data += '\n'
            report.write(data)

    def __init__(self, args, pipelinecommit, startingtime, scriptpath, analysistype, cutoff, pipeline, revbait):
        """
        :param args: command line arguments
        :param pipelinecommit: pipeline commit or version
        :param startingtime: time the script was started
        :param scriptpath: home path of the script
        :param analysistype: name of the analysis being performed - allows the program to find databases
        :param cutoff: percent identity cutoff for matches
        :param pipeline: boolean of whether this script needs to run as part of a particular assembly pipeline
        """
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.starttime = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        try:
            self.sequencepath = os.path.join(args.sequencepath, '')
        except AttributeError:
            self.sequencepath = self.path
        assert os.path.isdir(self.sequencepath), u'Sequence path  is not a valid directory {0!r:s}' \
            .format(self.sequencepath)
        try:
            self.targetpath = os.path.join(args.reffilepath, analysistype)
        except AttributeError:
            self.targetpath = os.path.join(args.targetpath, '')
        self.reportpath = os.path.join(self.path, 'reports')
        assert os.path.isdir(self.targetpath), u'Target path is not a valid directory {0!r:s}' \
            .format(self.targetpath)
        try:
            self.bcltofastq = args.bcltofastq
        except AttributeError:
            self.bcltofastq = False
        try:
            self.miseqpath = args.miseqpath
        except AttributeError:
            self.miseqpath = str()
        try:
            self.miseqfolder = args.miseqfolder
        except AttributeError:
            self.miseqfolder = str()
        try:
            self.fastqdestination = args.fastqdestination
        except AttributeError:
            self.fastqdestination = str()
        try:
            self.forwardlength = args.forwardlength
        except AttributeError:
            self.forwardlength = 'full'
        try:
            self.reverselength = args.reverselength
        except AttributeError:
            self.reverselength = 'full'
        self.numreads = 2 if self.reverselength != 0 else 1
        self.customsamplesheet = args.customsamplesheet
        self.logfile = args.logfile
        # Set the custom cutoff value
        self.cutoff = float(cutoff)
        try:
            self.averagedepth = int(args.averagedepth)
        except AttributeError:
            self.averagedepth = 10
        try:
            self.copy = args.copy
        except AttributeError:
            self.copy = False
        self.runmetadata = args.runmetadata
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        try:
            self.cpus = int(args.cpus)
        except AttributeError:
            self.cpus = multiprocessing.cpu_count()
        try:
            self.threads = int(self.cpus / len(self.runmetadata.samples)) if self.cpus / len(self.runmetadata.samples) \
                                                                             > 1 else 1
        except TypeError:
            self.threads = self.cpus
        self.taxonomy = {'Escherichia': 'coli', 'Listeria': 'monocytogenes', 'Salmonella': 'enterica'}
        self.analysistype = analysistype
        self.pipeline = pipeline
        self.revbait = revbait
        # Run the analyses
        self.runner()


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
    parser.add_argument('-d', '--fastqdestination',
                        help='Optional folder path to store .fastq files created using the fastqCreation module. '
                             'Defaults to path/miseqfolder')
    parser.add_argument('-r1', '--forwardlength',
                        default='full',
                        help='Length of forward reads to use. Can specify "full" to take the full length of '
                             'forward reads specified on the SampleSheet')
    parser.add_argument('-r2', '--reverselength',
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
    parser.add_argument('-a', '--averagedepth',
                        default=10,
                        help='Supply an integer of the minimum mapping depth in order to return a positive result ')
    parser.add_argument('-C', '--copy',
                        action='store_true',
                        help='Normally, the program will create symbolic links of the files into the sequence path, '
                             'however, the are occasions when it is necessary to copy the files instead')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.pipeline = False
    arguments.runmetadata.samples = MetadataObject()
    arguments.analysistype = 'genesippr'
    arguments.logfile = os.path.join(arguments.path, 'logfile')
    # Define the start time
    start = time.time()

    # Run the script
    GeneSippr(arguments, commit, start, homepath, arguments.analysistype, arguments.cutoff, arguments.pipeline, False)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')

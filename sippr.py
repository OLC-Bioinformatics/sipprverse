#!/usr/bin/python3
from accessoryFunctions.accessoryFunctions import make_path, MetadataObject, printtime
from accessoryFunctions.metadataprinter import MetadataPrinter
from spadespipeline.typingclasses import Resistance, Virulence
from sixteenS.sixteens_full import SixteenS as SixteensFull
from MLSTsippr.mlst import GeneSippr as MLSTSippr
from sipprverse_reporter.reports import Reports
from sipprCommon.objectprep import Objectprep
from sipprCommon.sippingmethods import Sippr
from serosippr.serosippr import SeroSippr
import MASHsippr.mash as mash
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
        if self.virulence:
            vir = Virulence(self, self.commit, self.starttime, self.homepath, 'virulence', 0.95, False, True)
            vir.reporter()
        if self.genesippr:
            # Run the genesippr analyses
            self.analysistype = 'genesippr'
            self.targetpath = os.path.join(self.reffilepath, self.analysistype, '')
            Sippr(self, 0.90)
            # Create the reports
            self.reports = Reports(self)
            Reports.reporter(self.reports)
        if self.sixteens:
            # Run the 16S analyses
            SixteensFull(self, self.commit, self.starttime, self.homepath, 'sixteens_full', 0.985)
        if self.mlst:
            '''
            Genus-specific
            '''
            MLSTSippr(self, self.commit, self.starttime, self.homepath, 'MLST', 1.0, True)
        if self.rmlst:
            MLSTSippr(self, self.commit, self.starttime, self.homepath, 'rMLST', 1.0, True)
        if self.resistance:
            # ResFinding
            res = Resistance(self, self.commit, self.starttime, self.homepath, 'resfinder', 0.9, False, True)
            res.main()
        if self.closestreference:
            self.pipeline = True
            mash.Mash(self, 'mash')
        if self.gdcs:
            # Run the GDCS analysis
            self.analysistype = 'GDCS'
            self.targetpath = os.path.join(self.targetpath, self.analysistype)
            Sippr(self, 0.95)
            # Create the reports
            Reports.gdcsreporter(self.reports)
        # Optionally perform serotyping
        if self.serotype:
            '''
            Genus-specific
            '''
            # # Perform serotyping for samples classified as Escherichia
            # for sample in self.runmetadata.samples:
            #     if sample.general.bestassemblyfile != 'NA':
            #         sample.mash = GenObject()
            #         try:
            #             sample.mash.closestrefseqgenus = sample.general.closestrefseqgenus
            #             for genus, species in self.taxonomy.items():
            #                 if genus == sample.mash.closestrefseqgenus:
            #                     sample.mash.closestrefseqspecies = species
            #         except KeyError:
            #             sample.mash.closestrefseqgenus = 'NA'
            #             sample.mash.closestrefseqspecies = 'NA'
            #     else:
            #         sample.mash.closestrefseqgenus = 'NA'
            #         sample.mash.closestrefseqspecies = 'NA'
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
        self.args = args
        # Define variables based on supplied arguments
        self.path = os.path.join(args.outputpath)
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.sequencepath = os.path.join(args.sequencepath)
        self.seqpath = self.sequencepath
        self.targetpath = os.path.join(args.referencefilepath)
        # ref file path is used to work with submodule code with a different naming scheme
        self.reffilepath = self.targetpath
        self.reportpath = os.path.join(self.path, 'reports')
        make_path(self.reportpath)
        assert os.path.isdir(self.targetpath), u'Target path is not a valid directory {0!r:s}' \
            .format(self.targetpath)
        # Set the custom cutoff value
        self.cutoff = args.customcutoffs
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.numthreads if args.numthreads else multiprocessing.cpu_count())
        self.closestreference = args.closestreference
        self.gdcs = args.gdcs
        self.genesippr = args.genesippr
        self.mlst = args.mlst
        self.resistance = args.resistance
        self.rmlst = args.rmlst
        self.serotype = args.serotype
        self.sixteens = args.sixteens
        self.virulence = args.virulence
        self.reports = str()
        self.threads = int()
        self.runmetadata = MetadataObject()
        self.taxonomy = {'Escherichia': 'coli', 'Listeria': 'monocytogenes', 'Salmonella': 'enterica'}
        self.analysistype = 'GeneSippr'
        self.pipeline = False
        self.logfile = os.path.join(self.path, 'log')


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
    parser.add_argument('-o', '--outputpath',
                        required=True,
                        help='Path to directory in which report folder is to be created')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    parser.add_argument('-r', '--referencefilepath',
                        help='Provide the location of the folder containing the pipeline accessory files (reference '
                             'genomes, MLST data, etc.')
    parser.add_argument('-n', '--numthreads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-u', '--customcutoffs',
                        default=0.90,
                        help='Custom cutoff values')
    parser.add_argument('-A', '--resistance',
                        action='store_true',
                        default=False,
                        help='Perform AMR analysis on samples')
    parser.add_argument('-C', '--closestreference',
                        action='store_true',
                        default=False,
                        help='Determine closest RefSeq match with mash')
    parser.add_argument('-G', '--genesippr',
                        action='store_true',
                        default=False,
                        help='Perform GeneSippr analysis on samples')
    parser.add_argument('-M', '--mlst',
                        action='store_true',
                        default=False,
                        help='Perform MLST analysis on samples')
    parser.add_argument('-Q', '--gdcs',
                        action='store_true',
                        default=False,
                        help='Perform GDCS Quality analysis on samples')
    parser.add_argument('-R', '--rmlst',
                        action='store_true',
                        default=False,
                        help='Perform rMLST analysis on samples')
    parser.add_argument('-S', '--serotype',
                        action='store_true',
                        default=False,
                        help='Perform serotype analysis on samples determined to be Escherichia')
    parser.add_argument('-V', '--virulence',
                        action='store_true',
                        default=False,
                        help='Perform virulence analysis on samples')
    parser.add_argument('-X', '--sixteens',
                        action='store_true',
                        default=False,
                        help='Perform 16S typing of samples')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    sippr = Sipprverse(arguments, commit, start, homepath)
    sippr.main()

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')

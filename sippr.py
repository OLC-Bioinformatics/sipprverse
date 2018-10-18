#!/usr/bin/python3
from accessoryFunctions.accessoryFunctions import make_path, MetadataObject, SetupLogging
from accessoryFunctions.metadataprinter import MetadataPrinter
from spadespipeline.typingclasses import Resistance, Virulence
from sixteenS.sixteens_full import SixteenS as SixteensFull
from MLSTsippr.mlst import GeneSippr as MLSTSippr
from customsippr.customsippr import CustomGenes
from sipprverse_reporter.reports import Reports
from sipprCommon.objectprep import Objectprep
from sipprCommon.sippingmethods import Sippr
from serosippr.serosippr import SeroSippr
import MASHsippr.mash as mash
from argparse import ArgumentParser
import multiprocessing
import subprocess
import logging
import time
import os

__author__ = 'adamkoziol'


class Sipprverse(object):

    def main(self):
        """
        Run the necessary methods in the correct order
        """
        logging.info('Starting {at} analysis pipeline'.format(at=self.analysistype))
        # Create the objects to be used in the analyses
        objects = Objectprep(self)
        objects.objectprep()
        self.runmetadata = objects.samples
        self.threads = int(self.cpus / len(self.runmetadata.samples)) if self.cpus / len(self.runmetadata.samples) > 1 \
            else 1
        if self.genesippr:
            # Run the genesippr analyses
            self.analysistype = 'genesippr'
            self.targetpath = os.path.join(self.reffilepath, self.analysistype)
            Sippr(self, 0.90)
            # Create the reports
            self.reports = Reports(self)
            Reports.reporter(self.reports)
        if self.sixteens:
            # Run the 16S analyses
            SixteensFull(args=self,
                         pipelinecommit=self.commit,
                         startingtime=self.starttime,
                         scriptpath=self.homepath,
                         analysistype='sixteens_full',
                         cutoff=0.985)
        if self.closestreference:
            self.pipeline = True
            mash.Mash(inputobject=self,
                      analysistype='mash')
        if self.rmlst:
            MLSTSippr(args=self,
                      pipelinecommit=self.commit,
                      startingtime=self.starttime,
                      scriptpath=self.homepath,
                      analysistype='rMLST',
                      cutoff=1.0,
                      pipeline=True)
        if self.resistance:
            # ResFinding
            res = Resistance(args=self,
                             pipelinecommit=self.commit,
                             startingtime=self.starttime,
                             scriptpath=self.homepath,
                             analysistype='resfinder',
                             cutoff=0.7,
                             pipeline=False,
                             revbait=True)
            res.main()
        if self.virulence:
            self.genus_specific()
            vir = Virulence(args=self,
                            pipelinecommit=self.commit,
                            startingtime=self.starttime,
                            scriptpath=self.homepath,
                            analysistype='virulence',
                            cutoff=0.95,
                            pipeline=False,
                            revbait=True)
            vir.reporter()
        if self.gdcs:
            # Run the GDCS analysis
            self.analysistype = 'GDCS'
            self.targetpath = os.path.join(self.reffilepath, self.analysistype)
            Sippr(inputobject=self,
                  cutoff=0.95)
            # Create the reports
            Reports.gdcsreporter(self.reports)
        if self.mlst:
            self.genus_specific()
            MLSTSippr(args=self,
                      pipelinecommit=self.commit,
                      startingtime=self.starttime,
                      scriptpath=self.homepath,
                      analysistype='MLST',
                      cutoff=1.0,
                      pipeline=True)
        # Optionally perform serotyping
        if self.serotype:
            self.genus_specific()
            SeroSippr(args=self,
                      pipelinecommit=self.commit,
                      startingtime=self.starttime,
                      scriptpath=self.homepath,
                      analysistype='serosippr',
                      cutoff=0.90,
                      pipeline=True)
        if self.user_genes:
            custom = CustomGenes(self)
            custom.main()
        # Print the metadata
        MetadataPrinter(self)

    def genus_specific(self):
        """
        For genus-specific targets, MLST and serotyping, determine if the closest refseq genus is known - i.e. if 16S
        analyses have been performed. Perform the analyses if required
        """
        # Initialise a variable to store whether the necessary analyses have already been performed
        closestrefseqgenus = False
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                try:
                    closestrefseqgenus = sample.general.closestrefseqgenus
                except AttributeError:
                    pass
        # Perform the 16S analyses as required
        if not closestrefseqgenus:
            logging.info('Must perform MASH analyses to determine genera of samples')
            # Run the analyses
            mash.Mash(self, 'mash')

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
        self.cpus = int(args.numthreads)
        self.closestreference = args.closestreference
        self.gdcs = args.gdcs
        self.genesippr = args.genesippr
        self.mlst = args.mlst
        self.resistance = args.resistance
        self.rmlst = args.rmlst
        self.serotype = args.serotype
        self.sixteens = args.sixteens
        self.virulence = args.virulence
        self.averagedepth = args.averagedepth
        try:
            self.user_genes = os.path.join(args.user_genes)
            assert os.path.isfile(self.user_genes), 'Cannot find user-supplied target file: {targets}. Please ' \
                                                    'double-check name and path of file'\
                .format(targets=self.user_genes)
        except TypeError:
            self.user_genes = args.user_genes
        # Set all the analyses to True if the full_suite option was selected
        if args.full_suite:
            self.closestreference = True
            self.gdcs = True
            self.genesippr = True
            self.mlst = True
            self.resistance = True
            self.rmlst = True
            self.serotype = True
            self.sixteens = True
            self.virulence = True
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
    parser = ArgumentParser(description='Performs GeneSipping on folder of FASTQ files')
    parser.add_argument('-o', '--outputpath',
                        required=True,
                        help='Path to directory in which report folder is to be created')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    parser.add_argument('-r', '--referencefilepath',
                        required=True,
                        help='Provide the location of the folder containing reference database')
    parser.add_argument('-a', '--averagedepth',
                        default=2,
                        help='Cutoff value for mapping depth to use when parsing BAM files.')
    parser.add_argument('-n', '--numthreads',
                        default=multiprocessing.cpu_count(),
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-c', '--customcutoffs',
                        default=0.90,
                        help='Custom cutoff values')
    parser.add_argument('-F', '--full_suite',
                        action='store_true',
                        default=False,
                        help='Perform all the built-in GeneSippr analyses (AMR, GDCS, Genesippr, MASH, MLST, '
                             'rMLST, Serotype, SixteenS, and Virulence')
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
    parser.add_argument('-U', '--user_genes',
                        default=False,
                        help='Name and path of user provided (multi-)FASTA file of genes to run against samples')
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
    SetupLogging()
    # Define the start time
    start = time.time()

    # Run the script
    sippr = Sipprverse(arguments, commit, start, homepath)
    sippr.main()

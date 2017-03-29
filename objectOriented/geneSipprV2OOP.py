#!/usr/bin/env python
# Import the necessary modules
# Subprocess->call is used for making system calls
import subprocess
import time
# Import the accessoryFunctions sub-module
import createFastq
import createObject
import SPAdesPipeline.OLCspades.metadataprinter as metadataprinter
from SPAdesPipeline.OLCspades.accessoryFunctions import *

__author__ = 'adamkoziol'


class GeneSippr(object):
    def objectprep(self):
        """
        Creates fastq files from an in-progress Illumina MiSeq run or create an object and moves files appropriately
        """
        printtime('Starting genesippr analysis pipeline', self.starttime)
        # Run the genesipping if necessary. Otherwise create the metadata object
        if self.bcltofastq:
            if self.customsamplesheet:
                assert os.path.isfile(self.customsamplesheet), 'Cannot find custom sample sheet as specified {}' \
                    .format(self.customsamplesheet)
            self.runmetadata = createFastq.FastqCreate(self)
        else:
            self.runmetadata = createObject.ObjectCreation(self)

    def runner(self):
        """
        Call the necessary methods in the appropriate order
        """
        import customtargets
        import sipprmash
        import sipprmlst
        import sixteenS
        # Create a sample object and create/link fastq files as necessary
        self.objectprep()
        if self.sixteens:
            sixteens = sixteenS.SixteenS(self, 'sixteens')
            sixteens.targets()
        metadataprinter.MetadataPrinter(self)
        # Run the typing modules
        if self.customtargetpath:
            custom = customtargets.Custom(self, 'custom', self.cutoff)
            custom.targets()
        else:
            #
            if self.rmlst:
                rmlst = sipprmlst.MLSTmap(self, 'rmlst')
                rmlst.targets()
            # Run the desired analyses
            sipprmash.SipprMash(self, 'mash')
            mlst = sipprmlst.MLSTmap(self, 'mlst')
            mlst.targets()

        metadataprinter.MetadataPrinter(self)

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: command line arguments
        :param pipelinecommit: pipeline commit or version
        :param startingtime: time the script was started
        :param scriptpath: home path of the script
        """
        import multiprocessing
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.starttime = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        self.args = args
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Output location is not a valid directory {0!r:s}'.format(self.path)
        self.sequencepath = os.path.join(args.sequencepath, '')
        assert os.path.isdir(self.sequencepath), u'Sequence path  is not a valid directory {0!r:s}' \
            .format(self.sequencepath)
        self.targetpath = os.path.join(args.targetpath, '')
        self.reportpath = os.path.join(self.path, 'reports')
        assert os.path.isdir(self.targetpath), u'Target path is not a valid directory {0!r:s}' \
            .format(self.targetpath)
        if args.customtargetpath:
            self.customtargetpath = os.path.join(args.customtargetpath, '')
            os.path.isdir(self.customtargetpath), u'Output location is not a valid directory {0!r:s}' \
                .format(self.customtargetpath)
        else:
            self.customtargetpath = ''
        self.bcltofastq = args.bcl2fastq
        self.fastqdestination = args.destinationfastq
        self.forwardlength = args.readlengthforward
        self.reverselength = args.readlengthreverse
        self.numreads = 2 if self.reverselength != 0 else 1
        self.customsamplesheet = args.customsamplesheet
        # Set the custom cutoff value
        self.cutoff = args.customcutoffs
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.numthreads if args.numthreads else multiprocessing.cpu_count())
        self.runmetadata = MetadataObject()
        #
        self.sixteens = args.sixteenStyping
        self.rmlst = args.rmlst
        # Run the analyses
        self.runner()


if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
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
    parser.add_argument('-16S', '--sixteenStyping',
                        action='store_true',
                        help='Perform 16S typing. Note that for analyses such as MLST, pathotyping, '
                             'serotyping, and virulence typing that require the genus of a strain to proceed, '
                             '16S typing will still be performed')
    parser.add_argument('-M', '--Mlst',
                        action='store_true',
                        help='Perform MLST analyses')
    parser.add_argument('-Y', '--pathotYping',
                        action='store_true',
                        help='Perform pathotyping analyses')
    parser.add_argument('-S', '--Serotyping',
                        action='store_true',
                        help='Perform serotyping analyses')
    parser.add_argument('-V',
                        '--Virulencetyping',
                        action='store_true',
                        help='Perform virulence typing analyses')
    parser.add_argument('-a', '--armi',
                        action='store_true',
                        help='Perform ARMI antimicrobial typing analyses')
    parser.add_argument('-r', '--rmlst',
                        action='store_true',
                        help='Perform rMLST analyses')
    parser.add_argument('-D', '--detailedReports',
                        action='store_true',
                        help='Provide detailed reports with percent identity and depth of coverage values '
                             'rather than just "+" for positive results')
    parser.add_argument('-C', '--customtargetpath',
                        help='Provide the path for a folder of custom targets .fasta format')
    parser.add_argument('-u', '--customcutoffs',
                        default=0.8,
                        help='Custom cutoff values')

    # TODO Add custom cutoffs
    # TODO Assert .fastq files present in provided folder
    # TODO Don't touch .fastq(.gz) files

    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    GeneSippr(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'

    # print json.dumps(seqdict, sort_keys=True, indent=4, separators=(',', ': '))
    """
    -m
    /media/miseq/MiSeqOutput
    -f
    160421_M02466_0152_000000000-AMLFW
    """
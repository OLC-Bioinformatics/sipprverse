#!/usr/bin/env python
import subprocess
import time
from sipprCommon.sippingmethods import *
from sipprCommon.objectprep import Objectprep
from accessoryFunctions.accessoryFunctions import *
from accessoryFunctions.metadataprinter import *

__author__ = 'adamkoziol'


class SixteenS(object):

    def runner(self):
        """
        Run the necessary methods in the correct order
        """
        printtime('Starting {} analysis pipeline'.format(self.analysistype), self.starttime)
        if not self.pipeline:
            # If the metadata has been passed from the method script, self.pipeline must still be false in order to
            # get Sippr() to function correctly, but the metadata shouldn't be recreated
            try:
                _ = vars(self.runmetadata)['samples']
            except KeyError:
                # Create the objects to be used in the analyses
                objects = Objectprep(self)
                objects.objectprep()
                self.runmetadata = objects.samples
            # Run the analyses
            Sippr(self, self.cutoff)
            #
            self.attributer()
            # Create the reports
            self.reporter()
            # Print the metadata
            printer = MetadataPrinter(self)
            printer.printmetadata()

    def attributer(self):
        """
        Parses the 16S target files to link accession numbers stored in the .fai and metadata files to the genera stored
        in the target file
        """
        from Bio import SeqIO
        import operator
        for sample in self.runmetadata.samples:
            # Load the records from the target file into a dictionary
            record_dict = SeqIO.to_dict(SeqIO.parse(sample[self.analysistype].baitfile, "fasta"))
            sample[self.analysistype].classification = set()
            sample[self.analysistype].genera = dict()
            # Add all the genera with hits into the set of genera
            for result in sample[self.analysistype].results:
                genus, species = record_dict[result].description.split('|')[-1].split()[:2]
                sample[self.analysistype].classification.add(genus)
                sample[self.analysistype].genera[result] = genus
            # Convert the set to a list for easier JSON serialisation
            sample[self.analysistype].classification = list(sample[self.analysistype].classification)
            # If there is a mixed sample, then further analyses will be complicated
            if len(sample[self.analysistype].classification) > 1:
                # print('multiple: ', sample.name, sample[self.analysistype].classification)
                sample.general.closestrefseqgenus = sample[self.analysistype].classification
                # sample.general.bestassemblyfile = 'NA'
                sample[self.analysistype].multiple = True
            else:
                sample[self.analysistype].multiple = False

                try:
                    # Recreate the results dictionary with the percent identity as a float rather than a string
                    sample[self.analysistype].intresults = \
                        {key: float(value) for key, value in sample[self.analysistype].results.items()}
                    # Set the best hit to be the top entry from the sorted results
                    sample[self.analysistype].besthit = sorted(sample[self.analysistype].intresults.items(),
                                                               key=operator.itemgetter(1), reverse=True)[0]
                    sample.general.closestrefseqgenus = sample[self.analysistype].classification[0]
                except IndexError:
                    sample.general.bestassemblyfile = 'NA'
                    # print('exception', sample.name, sample[self.analysistype].classification)

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        header = 'Strain,Gene,PercentIdentity,Genus,FoldCoverage\n'
        data = ''
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                if sample[self.analysistype].results:
                    if not sample[self.analysistype].multiple:
                        for name, identity in sample[self.analysistype].results.items():
                            if name == sample[self.analysistype].besthit[0]:
                                data += '{},{},{},{}\n'.format(name, identity, sample[self.analysistype].genera[name],
                                                               sample[self.analysistype].avgdepth[name])
                    else:
                        data += '{},{},{},{}\n'.format('multiple', 'NA', ';'.join(sample[self.analysistype]
                                                                                  .classification), 'NA')
                else:
                    data += '\n'
            report.write(header)
            report.write(data)

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
        self.analysistype = args.analysistype
        # Define variables based on supplied arguments
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.sequencepath = os.path.join(args.sequencepath, '')
        assert os.path.isdir(self.sequencepath), u'Sequence path  is not a valid directory {0!r:s}' \
            .format(self.sequencepath)
        self.targetpath = os.path.join(args.targetpath, self.analysistype, '')
        try:
            self.reportpath = args.reportpath
        except AttributeError:
            self.reportpath = os.path.join(self.path, 'reports')
        assert os.path.isdir(self.targetpath), u'Target path is not a valid directory {0!r:s}' \
            .format(self.targetpath)
        self.bcltofastq = args.bcltofastq
        self.miseqpath = args.miseqpath
        self.miseqfolder = args.miseqfolder
        self.fastqdestination = args.fastqdestination
        self.forwardlength = args.forwardlength
        self.reverselength = args.reverselength
        self.numreads = 2 if self.reverselength != 0 else 1
        self.customsamplesheet = args.customsamplesheet
        # Set the custom cutoff value
        self.cutoff = args.cutoff
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.cpus if args.cpus else multiprocessing.cpu_count())
        self.runmetadata = args.runmetadata
        self.taxonomy = {'Escherichia': 'coli', 'Listeria': 'monocytogenes', 'Salmonella': 'enterica'}
        self.pipeline = args.pipeline
        self.copy = args.copy
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
    parser.add_argument('-n', '--cpus',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-b', '--bcltofastq',
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
    parser.add_argument('-u', '--cutoff',
                        default=0.8,
                        help='Custom cutoff values')
    parser.add_argument('-C', '--copy',
                        action='store_true',
                        help='Normally, the program will create symbolic links of the files into the sequence path, '
                             'however, the are occasions when it is necessary to copy the files instead')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.pipeline = False
    arguments.runmetadata.samples = MetadataObject()
    arguments.analysistype = '16S'
    # Define the start time
    start = time.time()

    # Run the script
    SixteenS(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')

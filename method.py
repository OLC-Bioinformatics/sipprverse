#!/usr/bin/python3
from sipprCommon.sippingmethods import Sippr
from sipprCommon.objectprep import Objectprep
from accessoryFunctions.accessoryFunctions import MetadataObject, make_path, printtime
from accessoryFunctions.metadataprinter import MetadataPrinter
from sixteenS.sixteens_full import SixteenS as SixteensFull
from reporter.reports import Reports
from argparse import ArgumentParser
import multiprocessing
from time import sleep
from glob import glob
import subprocess
import time
import os
__author__ = 'adamkoziol'


class Method(object):

    def main(self):
        """
        Run the analyses using the inputted values for forward and reverse read length. However, if not all strains
        pass the quality thresholds, continue to periodically run the analyses on these incomplete strains until either
        all strains are complete, or the sequencing run is finished
        """
        printtime('Starting {} analysis pipeline'.format(self.analysistype), self.starttime, output=self.portallog)
        self.createobjects()
        # Run the genesipping analyses
        self.methods()
        # Determine if the analyses are complete
        self.complete()
        self.additionalsipping()
        # Update the report object
        self.reports = Reports(self)
        # Once all the analyses are complete, create reports for each sample
        Reports.methodreporter(self.reports)
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()

    def createobjects(self):
        # Set the name of the folders in which to store the current analysis based on the length of reads
        reads = '{}_{}'.format(self.forwardlength, self.reverselength)
        # Update the necessary variables to allow for the custom naming of folders based on the length forward and
        # reverse reads used to create the .fastq files
        self.fastqdestination = os.path.join(self.path, self.miseqfolder, reads)
        self.sequencepath = os.path.join(self.seqpath, reads)
        self.reportpath = os.path.join(self.reportpath, reads)
        self.samplesheetpath = os.path.join(self.path, 'SampleSheets', reads)
        # Create the objects to be used in the analyses
        objects = Objectprep(self)
        objects.objectprep()
        # Set the metadata
        self.runmetadata = objects.samples
        self.threads = int(self.cpus / len(self.runmetadata.samples)) if self.cpus / len(self.runmetadata.samples) > 1 \
            else 1
        # Pull the full length of the forward and reverse reads, as well as the indices
        self.forward = int(objects.forward)
        self.reverse = int(objects.reverse)
        self.index = objects.index
        # Store data from the sample sheet header and body
        self.header = objects.header
        self.rundata = objects.run
        # Set the list of all the sample names in the analysis
        self.samples = [sample.name for sample in self.runmetadata.samples]
        # If a custom sample sheet isn't specified, create one with all the samples
        if not self.customsamplesheet:
            for sample in self.runmetadata.samples:
                self.incomplete.append(sample.name)
            # Create the sample sheet
            self.samplesheet()
        # Set self.bcltofastq to False, as the call to Sippr() in self.methods will attempt to create the files again
        self.bcltofastq = False

    def additionalsipping(self):
        # If the analyses are not complete, continue to run the analyses until either all the strains pass the quality
        # thresholds, or until the sequencing run is complete
        while not self.analysescomplete:
            # Calculate the total number of reads needed for the run (forward + index1 + index2 + reverse). As the index
            # is the modified index used for the bcl2fastq its format is something like: AGGCAGAA-GCGTAAGA. Count only
            # the alphanumeric characters.
            self.sum = self.forward + sum(count.isalpha() for count in self.index) + self.reverse
            # Ensuring that the forward length is set to full - for testing only. In a real analysis, the forward length
            # has to be full due to the way that the sequencing is performed
            self.forwardlength = 'full'
            # Determine the number of completed cycles
            cycles = glob(
                os.path.join(self.miseqpath, self.miseqfolder, 'Data', 'Intensities', 'BaseCalls', 'L001', 'C*'))
            # If the run is complete, process the data one final time
            if len(cycles) == self.sum:
                printtime(
                    'Certain strains did not pass the quality thresholds. Final attempt of the pipeline. Using '
                    'the full reads'.format(self.forwardlength, self.reverselength, ),
                    self.starttime, output=self.portallog)
                # Set the boolean for the final iteration of the pipeline to true - this will allow all samples - even
                # ones considered incomplete to be entered into the final reports
                self.final = True
                # If the run is finished, then the reverse reads will be fully sequenced
                self.reverselength = 'full'
                # Set the name of the folders in which to store the current analysis based on the length of reads
                reads = '{}_{}'.format(self.forwardlength, self.reverselength)
                # Update the necessary variables to allow for the custom naming of folders based on the length forward
                # and reverse reads used to create the .fastq files
                self.fastqdestination = os.path.join(self.path, self.miseqfolder, reads)
                make_path(self.fastqdestination)
                self.sequencepath = os.path.join(self.seqpath, reads)
                make_path(self.sequencepath)
                self.reportpath = os.path.join(self.path, 'reports', reads)
                make_path(self.reportpath)
                self.samplesheetpath = os.path.join(self.path, 'SampleSheets', reads)
                # Create the sample sheet with only the samples that still need to be processed
                self.samplesheet()
                # Reset booleans
                self.analysescomplete = True
                self.bcltofastq = True
                # Create the objects to be used in the final analysis
                objects = Objectprep(self)
                objects.objectprep()
                # Set the metadata
                self.runmetadata = objects.samples
                self.bcltofastq = False
                # Run the analyses
                self.methods()
                self.complete()
            # If the sequencing run is not yet complete, continue to pull data from the MiSeq as it is created
            else:
                # Determine the length of reverse reads that can be used
                self.reverselength = str(len(cycles) - self.forward - sum(count.isalpha() for count in self.index))
                printtime(
                    'Certain strains did not pass the quality thresholds. Attempting the pipeline with the following '
                    'read lengths: forward {}, reverse {}'.format(self.forwardlength, self.reverselength),
                    self.starttime,
                    output=self.portallog)
                # Set the name of the folders in which to store the current analysis based on the length of reads
                reads = '{}_{}'.format(self.forwardlength, self.reverselength)
                # Update the necessary variables to allow for the custom naming of folders based on the length forward
                # and reverse reads used to create the .fastq files
                self.fastqdestination = os.path.join(self.path, self.miseqfolder, reads)
                make_path(self.fastqdestination)
                self.sequencepath = os.path.join(self.seqpath, reads)
                make_path(self.sequencepath)
                self.reportpath = os.path.join(self.path, 'reports', reads)
                make_path(self.reportpath)
                self.samplesheetpath = os.path.join(self.path, 'SampleSheets', reads)
                self.samplesheet()
                self.bcltofastq = True
                # Create the objects to be used in the analyses
                objects = Objectprep(self)
                objects.objectprep()
                # Set the metadata
                self.runmetadata = objects.samples
                self.methods()
                self.complete()
                # Allow the sequencer to complete approximately five cycles (~300 seconds per cycle) plus
                # however long it takes to run the analyses before trying again
                sleep(1500)

    def methods(self):
        self.run_genesippr()
        self.run_sixteens()
        self.run_gdcs()
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()

    def run_genesippr(self):
        # Run the genesippr analyses
        self.cutoff = 0.9
        self.analysistype = 'genesippr'
        self.targetpath = os.path.join(self.reffilepath, self.analysistype, '')
        Sippr(self, self.cutoff, 5)
        # Update the reports object
        self.reports = Reports(self)
        # Create the reports
        Reports.reporter(self.reports)
        Reports.genusspecific(self.reports)

    def run_sixteens(self):
        # Run the 16S analyses using the filtered database
        self.targetpath = self.reffilepath
        SixteensFull(self, self.commit, self.starttime, self.homepath, 'sixteens_full', 0.985)

    def run_gdcs(self):
        """

        """
        # Run the GDCS analysis
        self.analysistype = 'GDCS'
        self.pipeline = True
        Sippr(self, 0.95)
        # Create the reports
        Reports.gdcsreporter(self.reports)
        self.pipeline = False

    def complete(self):
        """
        Determine if the analyses of the strains are complete e.g. there are no missing GDCS genes, and the 
        sample.general.bestassemblyfile != 'NA'
        """
        # Boolean to store the completeness of the analyses
        allcomplete = True
        # Clear the list of samples that still require more sequence data
        self.incomplete = list()
        for sample in self.runmetadata.samples:
            try:
                # If the sample has been tagged as incomplete, only add it to the complete metadata list if the
                # pipeline is on its final iteration
                if sample.general.incomplete:
                    if self.final:
                        self.completemetadata.append(sample)
                    else:
                        sample.general.complete = False
                        allcomplete = False
                        self.incomplete.append(sample.name)
            except KeyError:
                sample.general.complete = True
                self.completemetadata.append(sample)
        # If all the samples are complete, set the global variable for run completeness to True
        if allcomplete:
            self.analysescomplete = True

    def samplesheet(self):
        """
        Create a custom sample sheet based on the original sample sheet for the run, but only including the samples
        that did not pass the quality threshold on the previous iteration
        """
        make_path(self.samplesheetpath)
        self.customsamplesheet = os.path.join(self.samplesheetpath, 'SampleSheet.csv')
        header = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID',
                  'index2', 'Sample_Project', 'Description']
        with open(self.customsamplesheet, 'w') as samplesheet:
            lines = str()
            lines += '[Header]\n'
            lines += 'IEMFileVersion,{}\n'.format(self.header['IEMFileVersion'])
            lines += 'Investigator Name,{}\n'.format(self.header['InvestigatorName'])
            lines += 'Experiment Name,{}\n'.format(self.header['ExperimentName'])
            lines += 'Date,{}\n'.format(self.header['Date'])
            lines += 'Workflow,{}\n'.format(self.header['Workflow'])
            lines += 'Application,{}\n'.format(self.header['Application'])
            lines += 'Assay,{}\n'.format(self.header['Assay'])
            lines += 'Description,{}\n'.format(self.header['Description'])
            lines += 'Chemistry,{}\n'.format(self.header['Chemistry'])
            lines += '\n'
            lines += '[Reads]\n'
            lines += str(self.forward) + '\n'
            lines += str(self.reverse) + '\n'
            lines += '\n'
            lines += '[Settings]\n'
            lines += 'ReverseComplement,{}\n'.format(self.header['ReverseComplement'])
            lines += 'Adapter,{}\n'.format(self.header['Adapter'])
            lines += '\n'
            lines += '[Data]\n'
            lines += ','.join(header)
            lines += '\n'
            # Correlate all the samples added to the list of incomplete samples with their metadata
            for incomplete in self.incomplete:
                for sample in self.rundata:
                    if incomplete == sample['SampleID']:
                        # Use each entry in the header list as a key for the rundata dictionary
                        for data in header:
                            # Modify the key to be consistent with how the dictionary was populated
                            result = sample[data.replace('_', '')]
                            # Description is the final entry in the list, and shouldn't have a , following the value
                            if data != 'Description':
                                lines += '{},'.format(result.replace('NA', ''))
                            # This entry should have a newline instead of a ,
                            else:
                                lines += '{}\n'.format(result.replace('NA', ''))
            # Write the string to the sample sheet
            samplesheet.write(lines)

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
        self.path = os.path.join(args.outputpath)
        make_path(self.path)
        assert os.path.isdir(self.path), 'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        try:
            self.portallog = args.portallog
        except AttributeError:
            self.portallog = os.path.join(self.path, 'portal.log')
        try:
            os.remove(self.portallog)
        except FileNotFoundError:
            pass
        self.sequencepath = os.path.join(self.path, 'sequences')
        self.seqpath = self.sequencepath
        self.targetpath = os.path.join(args.targetpath)
        # ref file path is used to work with sub module code with a different naming scheme
        self.reffilepath = self.targetpath
        self.reportpath = os.path.join(self.path, 'reports')
        make_path(self.reportpath)
        assert os.path.isdir(self.targetpath), 'Target path is not a valid directory {0!r:s}' \
            .format(self.targetpath)
        self.bcltofastq = True
        self.miseqpath = args.miseqpath
        self.miseqfolder = args.miseqfolder
        # self.fastqdestination = args.destinationfastq
        self.fastqdestination = str()
        self.forwardlength = args.readlengthforward
        self.reverselength = args.readlengthreverse
        self.numreads = 2 if self.reverselength != 0 else 1
        self.customsamplesheet = args.customsamplesheet
        # Set the custom cutoff value
        self.cutoff = float()
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        try:
            self.cpus = int(args.numthreads)
        except (AttributeError, TypeError):
            self.cpus = multiprocessing.cpu_count()
        self.threads = int()
        self.runmetadata = MetadataObject()
        self.taxonomy = {'Escherichia': 'coli', 'Listeria': 'monocytogenes', 'Salmonella': 'enterica'}
        self.analysistype = 'GeneSipprMethod'
        self.copy = args.copy
        try:
            self.debug = args.debug
        except AttributeError:
            self.debug = False
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
    parser.add_argument('-r', '--referencefilepath',
                        help='Provide the location of the folder containing the pipeline accessory files (reference '
                             'genomes, MLST data, etc.')
    parser.add_argument('-n', '--numthreads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-b', '--bcl2fastq',
                        action='store_true',
                        help='Optionally run bcl2fastq on an in-progress Illumina MiSeq run. Must include:'
                             'miseqpath, and miseqfolder arguments, and optionally readlengthforward, '
                             'readlengthreverse, and projectName arguments.')
    parser.add_argument('-m', '--miseqpath',
                        required=True,
                        help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', '--miseqfolder',
                        required=True,
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
    parser.add_argument('-C', '--copy',
                        action='store_true',
                        help='Normally, the program will create symbolic links of the files into the sequence path, '
                             'however, the are occasions when it is necessary to copy the files instead')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.portallog = os.path.join(arguments.path, 'portal.log')
    # Define the start time
    start = time.time()

    # Run the script
    method = Method(arguments, commit, start, homepath)
    method.main()

    # Print a bold, green exit statement
    printtime('Analyses complete', start, option='\033[1;92m', output=arguments.portallog)

"""
-b
-m
/media/miseq
-f
170420_M02466_0033_000000000-AYMM9
-r1
120
-r2
0
-C


/home/adamkoziol/Bioinformatics/sippr/method
-t
/mnt/nas/bio_requests/8312/newsixteens/targets
-b
-m
/home/adamkoziol/Bioinformatics/
-f
161104_M02466_0002_000000000-AV4G5
-r1
20
-r2
0
-C
-c
/home/adamkoziol/Bioinformatics/sippr/method/SampleSheet.csv


-b
-m
/home/adamkoziol/Bioinformatics/
-f
161104_M02466_0002_000000000-AV4G5
-r1
full
-r2
full
-c
/home/adamkoziol/Bioinformatics/sippr/method/SampleSheet.csv

"""

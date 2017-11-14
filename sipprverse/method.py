#!/usr/bin/python3
from sipprCommon.sippingmethods import Sippr
from sipprCommon.objectprep import Objectprep
from accessoryFunctions.accessoryFunctions import MetadataObject, make_path, printtime
from accessoryFunctions.metadataprinter import *
from sixteenS.sixteens_full import SixteenS as SixteensFull
from argparse import ArgumentParser
from glob import glob
import subprocess
import numpy
import time
import os
__author__ = 'adamkoziol'


class Method(object):

    def method(self):
        """
        Run the analyses using the inputted values for forward and reverse read length. However, if not all strains
        pass the quality thresholds, continue to periodically run the analyses on these incomplete strains until either
        all strains are complete, or the sequencing run is finished
        """
        printtime('Starting {} analysis pipeline'.format(self.analysistype), self.starttime)
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
        # Run the genesipping analyses
        self.methods()
        # Determine if the analyses are complete
        self.complete()
        # Calculate the total number of reads required for the run (forward + index1 + index2 + reverse). As the index
        # is the modified index used for the bcl2fastq its format is something like: AGGCAGAA-GCGTAAGA. Count only
        # the alphanumeric characters.
        self.sum = self.forward + sum(count.isalpha() for count in self.index) + self.reverse
        # If the analyses are not complete, continue to run the analyses until either all the strains pass the quality
        # thresholds, or until the sequencing run is complete
        while not self.analysescomplete:
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
                    'the full reads'.format(self.forwardlength, self.reverselength),
                    self.starttime)
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
                make_path(self.seqpath)
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
                from time import sleep
                printtime(
                    'Certain strains did not pass the quality thresholds. Attempting the pipeline with the following'
                    'read lengths: forward {}, reverse {}'.format(self.forwardlength, self.reverselength),
                    self.starttime)
                # Determine the length of reverse reads that can be used
                self.reverselength = self.sum - len(cycles)
                # Set the name of the folders in which to store the current analysis based on the length of reads
                reads = '{}_{}'.format(self.forwardlength, self.reverselength)
                # Update the necessary variables to allow for the custom naming of folders based on the length forward
                # and reverse reads used to create the .fastq files
                self.fastqdestination = os.path.join(self.path, self.miseqfolder, reads)
                make_path(self.fastqdestination)
                self.sequencepath = os.path.join(self.seqpath, reads)
                make_path(self.seqpath)
                self.reportpath = os.path.join(self.path, 'reports', reads)
                make_path(self.reportpath)
                self.samplesheetpath = os.path.join(self.path, 'SampleSheets', reads)
                self.samplesheet()
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
        # Once all the analyses are complete, create reports for each sample
        self.methodreporter()
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()

    def methods(self):
        """
        Method to allow the analyses to be called in a repeatable fashion 
        """
        # Run the genesippr analyses
        self.cutoff = 0.8
        self.analysistype = 'genesippr'
        self.targetpath = os.path.join(self.reffilepath, self.analysistype, '')
        Sippr(self, self.cutoff)
        # Create the reports
        self.reporter()
        # Run the 16S analyses using the filtered database
        self.targetpath = self.reffilepath
        SixteensFull(self, self.commit, self.starttime, self.homepath, 'sixteens_full', 0.985)
        # Run the GDCS analysis
        self.analysistype = 'GDCS'
        self.pipeline = True
        Sippr(self, 0.95)
        # Create the reports
        self.gdcsreporter()
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()
        self.pipeline = False

    def runner(self):
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
        self.reporter()
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
        self.gdcsreporter()
        '''
        from serosippr.serosippr import SeroSippr
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
        '''
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()

    def reporter(self):
        """
        Creates a report of the genesippr results
        """
        from Bio import SeqIO
        printtime('Creating {} report'.format(self.analysistype), self.starttime)
        # Create a dictionary to link all the genera with their genes
        genusgenes = dict()
        # A list to store all the unique gene names
        geneset = list()
        # The organism-specific targets are in .tfa files in the target path
        targetpath = str()
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                targetpath = sample[self.analysistype].targetpath
        for organismfile in glob(os.path.join(targetpath, '*.tfa')):
            organism = os.path.splitext(os.path.basename(organismfile))[0]
            # Use BioPython to extract all the gene names from the file
            for record in SeqIO.parse(open(organismfile), 'fasta'):
                # Add the gene name to the list of genes if it is not already present. I wanted to use a set, but
                # I also wanted to keep the input order, which is why I used the if .. not in loop instead
                if record.id not in geneset:
                    geneset.append(record.id)
                # Append the gene names to the genus-specific list
                try:
                    genusgenes[organism].append(record.id)
                except (KeyError, IndexError):
                    genusgenes[organism] = list()
                    genusgenes[organism].append(record.id)
        # Determine from which genera the gene hits were sourced
        for sample in self.runmetadata.samples:
            # Initialise the list to store the genera
            sample[self.analysistype].targetgenera = list()
            if sample.general.bestassemblyfile != 'NA':
                for organism in genusgenes:
                    # Iterate through all the genesippr hits and attribute each gene to the appropriate genus
                    for gene in sample[self.analysistype].results:
                        # If the gene name is in the genes from that organism, add the genus name to the list of
                        # genera found in the sample
                        if gene in genusgenes[organism]:
                            if organism not in sample[self.analysistype].targetgenera:
                                sample[self.analysistype].targetgenera.append(organism)
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # The report will have every gene for all genera in the header
        header = 'Strain,Genus,{},\n'.format(','.join(geneset))
        data = str()
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                if sample.general.bestassemblyfile != 'NA':
                    # Add the genus/genera found in the sample
                    data += '{},{},'.format(sample.name, ';'.join(sample[self.analysistype].targetgenera))
                    if sample[self.analysistype].results:
                        for gene in geneset:
                            # If the gene was not found in the sample, print an empty cell in the report
                            if gene not in [target[0] for target in sample[self.analysistype].results.items()]:
                                data += ','
                            # Print the required information for the gene
                            for name, identity in sample[self.analysistype].results.items():
                                if name == gene:
                                    data += '{}% ({} +/- {}),'.format(identity,
                                                                      sample[self.analysistype].avgdepth[name],
                                                                      sample[self.analysistype].standarddev[name])
                        # Add a newline after each sample
                        data += '\n'
                    # Add a newline if the sample did not have any gene hits
                    else:
                        data += '\n'
            # Write the header and data to file
            report.write(header)
            report.write(data)

    def gdcsreporter(self):
        """
        Creates a report of the GDCS results
        """
        printtime('Creating {} report'.format(self.analysistype), self.starttime)
        # Initialise list to store all the GDCS genes, and genera in the analysis
        gdcs = list()
        genera = list()
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                sample[self.analysistype].createreport = True
                # Determine which genera are present in the analysis
                if sample.general.closestrefseqgenus not in genera:
                    genera.append(sample.general.closestrefseqgenus)
                try:
                    # Add all the GDCS genes to the list
                    for gene in sorted(sample[self.analysistype].faidict):
                        if gene not in gdcs:
                            gdcs.append(gene)
                except KeyError:
                    sample[self.analysistype].createreport = False
            else:
                sample[self.analysistype].createreport = False
        header = 'Strain,Genus,Matches,MeanCoverage,Pass/Fail,{},\n'.format(','.join(gdcs))
        data = str()
        with open(os.path.join(self.reportpath, '{}.csv'.format(self.analysistype)), 'w') as report:
            # Sort the samples in the report based on the closest refseq genus e.g. all samples with the same genus
            # will be grouped together in the report
            for genus in genera:
                for sample in self.runmetadata.samples:
                    if sample.general.closestrefseqgenus == genus:
                        if sample[self.analysistype].createreport:
                            sample[self.analysistype].totaldepth = list()
                            # Add the sample to the report if it matches the current genus
                            # if genus == sample.general.closestrefseqgenus:
                            data += '{},{},'.format(sample.name, genus)
                            # Initialise a variable to store the number of GDCS genes were matched
                            count = 0
                            # As I want the count to be in the report before all the gene results, this string will
                            # store the specific sample information, and will be added to data once count is known
                            specific = str()
                            for gene in gdcs:
                                # As there are different genes present in the GDCS databases for each organism of
                                # interest, genes that did not match because they're absent in the specific database are
                                # indicated using an X
                                if gene not in [result for result in sample[self.analysistype].faidict]:
                                    specific += 'X,'
                                else:
                                    try:
                                        # Report the necessary information for each gene result
                                        identity = sample[self.analysistype].results[gene]
                                        specific += '{}% ({} +/- {}),'\
                                            .format(identity, sample[self.analysistype].avgdepth[gene],
                                                    sample[self.analysistype].standarddev[gene])
                                        sample[self.analysistype].totaldepth.append(
                                            float(sample[self.analysistype].avgdepth[gene]))
                                        count += 1
                                    # If the gene was missing from the results attribute, add a - to the cell
                                    except KeyError:
                                        sample.general.incomplete = True
                                        specific += '-,'
                            # Calculate the mean depth of the genes and the standard deviation
                            sample[self.analysistype].mean = numpy.mean(sample[self.analysistype].totaldepth)
                            sample[self.analysistype].stddev = numpy.std(sample[self.analysistype].totaldepth)
                            # Determine whether the sample pass the necessary quality criteria:
                            # Pass, all GDCS, mean coverage greater than 20X coverage;
                            # ?: Indeterminate value;
                            # -: Fail value
                            if count == len(sample[self.analysistype].faidict):
                                if sample[self.analysistype].mean > 20:
                                    quality = '+'
                                else:
                                    quality = '?'
                            else:
                                quality = '-'
                            # Add the count, mean depth with standard deviation, the pass/fail determination,
                            #  and the total number of GDCS genes as well as the results
                            data += '{hits}/{total},{mean} +/- {std},{fail},{gdcs}\n'\
                                .format(hits=str(count),
                                        total=len(sample[self.analysistype].faidict),
                                        mean='{:.2f}'.format(sample[self.analysistype].mean),
                                        std='{:.2f}'.format(sample[self.analysistype].stddev),
                                        fail=quality,
                                        gdcs=specific)
                        # Any samples with a best assembly of 'NA' are considered incomplete.
                        else:
                            data += '{},{},,,-\n'.format(sample.name, sample.general.closestrefseqgenus)
                            sample.general.incomplete = True
            # Write the header and data to file
            report.write(header)
            report.write(data)

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

    def methodreporter(self):
        """
        Create final reports collating results from all the individual iterations through the method pipeline
        """
        # Ensure that the analyses are set to complete
        self.analysescomplete = True
        # Reset the report path to original value
        self.reportpath = os.path.join(self.path, 'reports')
        # Clear the runmetadata - it will be populated with all the metadata from completemetadata
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = list()
        # As the samples were entered into self.completemetadata depending on when they passed the quality threshold,
        # this list is not ordered numerically/alphabetically like the original runmetadata. Reset the order.
        for strain in self.samples:
            for sample in self.completemetadata:
                if sample.name == strain:
                    # Append the sample to the ordered list of objects
                    self.runmetadata.samples.append(sample)
        # Set the analysis type for each analysis performed
        self.analysistype = 'genesippr'
        # Create the reports
        self.reporter()
        self.analysistype = '16S'
        self.sixteensreporter()
        # Run the GDCS analysis
        self.analysistype = 'GDCS'
        # Create the reports
        self.gdcsreporter()

    def sixteensreporter(self):
        """
        Creates a report of the results
        """
        printtime('Creating {} report'.format(self.analysistype), self.starttime)
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
                                data += '{},{},{},{}\n'.format(name, identity,
                                                               sample[self.analysistype].genera[name],
                                                               sample[self.analysistype].avgdepth[name])
                    else:
                        data += '{},{},{},{}\n'.format('multiple', 'NA',
                                                       ';'.join(sample[self.analysistype].classification), 'NA')
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
        self.cutoff = float(args.customcutoffs)
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.numthreads if args.numthreads else multiprocessing.cpu_count())
        self.threads = int()
        self.runmetadata = MetadataObject()
        self.taxonomy = {'Escherichia': 'coli', 'Listeria': 'monocytogenes', 'Salmonella': 'enterica'}
        self.analysistype = 'GeneSippr method'
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
        if self.bcltofastq:
            make_path(self.sequencepath)
            self.method()
        else:
            # Run the analyses
            assert os.path.isdir(self.sequencepath), u'Sequence path  is not a valid directory {0!r:s}' \
                .format(self.sequencepath)
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
    Method(arguments, commit, start, homepath)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')

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

"""

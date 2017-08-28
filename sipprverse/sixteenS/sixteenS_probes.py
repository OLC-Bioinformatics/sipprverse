#!/usr/bin/env python 3

import subprocess
import time
from sipprCommon.sippingmethods import *
from sipprCommon.objectprep import Objectprep
from accessoryFunctions.accessoryFunctions import *
from accessoryFunctions.metadataprinter import *

__author__ = 'adamkoziol'


class ProbeSippr(Sippr):

    def targets(self):
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].targetpath = self.targetpath
                baitpath = os.path.join(self.targetpath, 'bait')
                sample[self.analysistype].baitfile = glob(os.path.join(baitpath, '*.fa'))[0]
                # Create the hash file of the baitfile
                targetbase = sample[self.analysistype].baitfile.split('.')[0]
                sample[self.analysistype].hashfile = targetbase + '.mhs.gz'
                sample[self.analysistype].hashcall = 'cd {} && mirabait -b {} -k 19 -K {}' \
                    .format(sample[self.analysistype].targetpath,
                            sample[self.analysistype].baitfile,
                            sample[self.analysistype].hashfile)
                if not os.path.isfile(sample[self.analysistype].hashfile):
                    call(sample[self.analysistype].hashcall, shell=True, stdout=self.devnull, stderr=self.devnull)
                # Ensure that the hash file was successfully created
                assert os.path.isfile(sample[self.analysistype].hashfile), \
                    u'Hashfile could not be created for the target file {0!r:s}'.format(
                        sample[self.analysistype].baitfile)
                sample[self.analysistype].outputdir = os.path.join(sample.run.outputdirectory, self.analysistype)
                sample[self.analysistype].baitedfastq = \
                    '{}/{}_targetMatches.fastq'.format(sample[self.analysistype].outputdir, self.analysistype)
                sample[self.analysistype].phylogeny = list()
                sample[self.analysistype].complete = False
        #
        self.baiting()

    def baiting(self):
        # Perform baiting
        printtime('Performing kmer baiting of fastq files with {} targets'.format(self.analysistype), self.start)
        # Create and start threads for each fasta file in the list
        for i in range(len(self.runmetadata)):
            # Send the threads to the bait method
            threads = Thread(target=self.bait, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                # Add the sample to the queue
                self.baitqueue.put(sample)
        self.baitqueue.join()
        #
        self.premap()

    def premap(self):
        complete = False
        incomplete = list()
        analysis = str()
        while not complete:
            for sample in self.runmetadata:
                if sample.general.bestassemblyfile != 'NA':
                    if not sample[self.analysistype].complete:
                        try:
                            #
                            currentpath = os.path.join(sample[self.analysistype].targetpath,
                                                       *sample[self.analysistype].phylogeny)

                            currentpath = os.path.join(currentpath, '')
                            currenttarget = glob(currentpath + '*.fa')[0]
                            base = os.path.basename(currenttarget).split('_')[0]
                            currentanalysis = '{}_{}'.format(self.analysistype, base)
                            analysis = currentanalysis
                            setattr(sample, currentanalysis, GenObject())
                            #
                            sample[currentanalysis].outputdir = os.path.join(sample[self.analysistype].outputdir, base)
                            make_path(sample[currentanalysis].outputdir)
                            sample[currentanalysis].baitfile = currenttarget
                            sample[currentanalysis].baitedfastq = sample[self.analysistype].baitedfastq
                            sample[currentanalysis].hashfile = sample[self.analysistype].hashfile
                            sample[currentanalysis].targetpath = currentpath
                            incomplete.append(sample)
                        except IndexError:
                            print('error', sample.name, sample[self.analysistype].phylogeny)
                            pass
            if incomplete:
                self.mapping(analysis, incomplete)
                self.indexing(analysis, incomplete)
                self.parsing(analysis, incomplete)
                self.postmapping(analysis, incomplete)
                # complete = True
            else:
                complete = True
            incomplete = list()
            # print(sample[self.analysistype].datastore)
        self.reporting()

    def mapping(self, analysistype, metadata):
        """

        :param analysistype:
        :param metadata:
        """
        printtime('Performing reference mapping', self.start)
        for i in range(len(metadata)):
            # Send the threads to
            threads = Thread(target=self.map, args=())
            # Set the daemon to True - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Set the path/name for the sorted bam file to be created
                sample[analysistype].sortedbam = '{}/{}_sorted.bam'.format(sample[analysistype].outputdir,
                                                                           analysistype)
                # Remove the file extension of the bait file for use in the indexing command
                sample[analysistype].baitfilenoext = sample[analysistype].baitfile.split('.')[0]
                # Use bowtie2 wrapper to create index the target file
                bowtie2build = Bowtie2BuildCommandLine(reference=sample[analysistype].baitfile,
                                                       bt2=sample[analysistype].baitfilenoext,
                                                       **self.builddict)
                # Use samtools wrapper to set up the bam sorting command
                ''''''
                samsort = SamtoolsSortCommandline(input=sample[analysistype].sortedbam,
                                                  o=True,
                                                  out_prefix="-")
                samtools = [
                    # When bowtie2 maps reads to all possible locations rather than choosing a 'best' placement, the
                    # SAM header for that read is set to 'secondary alignment', or 256. Please see:
                    # http://davetang.org/muse/2014/03/06/understanding-bam-flags/ The script below reads in the stdin
                    # and subtracts 256 from headers which include 256
                    'python {}/sipprcommon/editsamheaders.py'.format(self.homepath),
                    # Use samtools wrapper to set up the samtools view
                    SamtoolsViewCommandline(b=True,
                                            S=True,
                                            h=True,
                                            input_file="-"),
                    samsort]
                # Add custom parameters to a dictionary to be used in the bowtie2 alignment wrapper
                indict = {'--very-sensitive-local': True,
                          # For short targets, the match bonus can be increased
                          '--ma': self.matchbonus,
                          '-U': sample[analysistype].baitedfastq,
                          '-a': True,
                          '--threads': self.cpus,
                          '--local': True}
                # Create the bowtie2 reference mapping command
                bowtie2align = Bowtie2CommandLine(bt2=sample[analysistype].baitfilenoext,
                                                  threads=self.cpus,
                                                  samtools=samtools,
                                                  **indict)
                # Create the command to faidx index the bait file
                sample[analysistype].faifile = sample[analysistype].baitfile + '.fai'
                samindex = SamtoolsFaidxCommandline(reference=sample[analysistype].baitfile)
                # Add the commands (as strings) to the metadata
                sample[analysistype].samindex = str(samindex)
                # Add the commands to the queue. Note that the commands would usually be set as attributes of the sample
                # but there was an issue with their serialization when printing out the metadata
                if not os.path.isfile(sample[analysistype].baitfilenoext + '.1' + self.bowtiebuildextension):
                    stdoutbowtieindex, stderrbowtieindex = map(StringIO,
                                                               bowtie2build(cwd=sample[analysistype].targetpath))
                    # Write any error to a log file
                    if stderrbowtieindex:
                        # Write the standard error to log, bowtie2 puts alignment summary here
                        with open(os.path.join(sample[analysistype].targetpath,
                                               '{}_bowtie_index.log'.format(analysistype)), 'a+') as log:
                            log.writelines(logstr(bowtie2build, stderrbowtieindex.getvalue(),
                                                  stdoutbowtieindex.getvalue()))
                    # Close the stdout and stderr streams
                    stdoutbowtieindex.close()
                    stderrbowtieindex.close()
                self.mapqueue.put((sample, bowtie2build, bowtie2align, samindex, analysistype))
        self.mapqueue.join()

    def map(self):
        while True:
            # Get the necessary values from the queue
            sample, bowtie2build, bowtie2align, samindex, analysistype = self.mapqueue.get()
            # Use samtools faidx to index the bait file - this will be used in the sample parsing
            if not os.path.isfile(sample[analysistype].faifile):
                stdoutindex, stderrindex = map(StringIO, samindex(cwd=sample[analysistype].targetpath))
                # Write any error to a log file
                if stderrindex:
                    # Write the standard error to log, bowtie2 puts alignment summary here
                    with open(os.path.join(sample[analysistype].targetpath,
                                           '{}_samtools_index.log'.format(analysistype)), 'a+') as log:
                        log.writelines(logstr(samindex, stderrindex.getvalue(), stdoutindex.getvalue()))
                # Close the stdout and stderr streams
                stdoutindex.close()
                stderrindex.close()
            # Only run the functions if the sorted bam files and the indexed bait file do not exist
            if not os.path.isfile(sample[analysistype].sortedbam):
                # Set stdout to a stringIO stream
                stdout, stderr = map(StringIO, bowtie2align(cwd=sample[analysistype].outputdir))
                if stderr:
                    # Write the standard error to log, bowtie2 puts alignment summary here
                    with open(os.path.join(sample[analysistype].outputdir,
                                           '{}_bowtie_samtools.log'.format(analysistype)), 'a+') as log:
                        log.writelines(logstr([bowtie2align], stderr.getvalue(), stdout.getvalue()))
                stdout.close()
                stderr.close()
            self.mapqueue.task_done()

    def indexing(self, analysistype, metadata):
        printtime('Indexing sorted bam files', self.start)
        for i in range(len(metadata)):
            # Send the threads to
            threads = Thread(target=self.index, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in metadata:
            if sample.general.bestassemblyfile != 'NA':
                bamindex = SamtoolsIndexCommandline(input=sample[analysistype].sortedbam)
                sample[analysistype].sortedbai = sample[analysistype].sortedbam + '.bai'
                sample[analysistype].bamindex = str(bamindex)
                self.indexqueue.put((sample, bamindex, analysistype))
        self.indexqueue.join()

    def index(self):
        while True:
            sample, bamindex, analysistype = self.indexqueue.get()
            # Only make the call if the .bai file doesn't already exist
            if not os.path.isfile(sample[analysistype].sortedbai):
                # Use cStringIO streams to handle bowtie output
                stdout, stderr = map(StringIO, bamindex(cwd=sample[analysistype].outputdir))
                if stderr:
                    # Write the standard error to log
                    with open(os.path.join(sample[analysistype].outputdir,
                                           '{}_samtools_bam_index.log'.format(analysistype)), 'a+') as log:
                        log.writelines(logstr(bamindex, stderr.getvalue(), stdout.getvalue()))
                stderr.close()
            self.indexqueue.task_done()

    def parsing(self, analysistype, metadata):
        printtime('Parsing sorted bam files', self.start)
        for i in range(len(metadata)):
            # Send the threads to
            threads = Thread(target=self.parse, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Get the fai file into a dictionary to be used in parsing results
                with open(sample[analysistype].faifile, 'r') as faifile:
                    for line in faifile:
                        data = line.split('\t')
                        try:
                            sample[analysistype].faidict[data[0]] = int(data[1])
                        except KeyError:
                            sample[analysistype].faidict = dict()
                            sample[analysistype].faidict[data[0]] = int(data[1])
                self.parsequeue.put((sample, analysistype))
        self.parsequeue.join()

    def parse(self):
        import pysamstats
        import operator
        import numpy
        while True:
            sample, analysistype = self.parsequeue.get()
            # Initialise dictionaries to store parsed data
            matchdict = dict()
            depthdict = dict()
            seqdict = dict()
            snpdict = dict()
            gapdict = dict()
            maxdict = dict()
            mindict = dict()
            deviationdict = dict()
            sample[analysistype].results = dict()
            sample[analysistype].avgdepth = dict()
            sample[analysistype].resultssnp = dict()
            sample[analysistype].resultsgap = dict()
            sample[analysistype].sequences = dict()
            sample[analysistype].maxcoverage = dict()
            sample[analysistype].mincoverage = dict()
            sample[analysistype].standarddev = dict()
            # Variable to store the expected position in gene/allele
            pos = 0
            try:
                # Use the stat_variation function of pysam stats to return records parsed from sorted bam files
                # Values of interest can be retrieved using the appropriate keys
                for rec in pysamstats.stat_variation(alignmentfile=sample[analysistype].sortedbam,
                                                     fafile=sample[analysistype].baitfile,
                                                     max_depth=1000000):
                    # Initialise seqdict with the current gene/allele if necessary with an empty string
                    if rec['chrom'] not in seqdict:
                        seqdict[rec['chrom']] = str()
                        # Since this is the first position in a "new" gene/allele, reset the pos variable to 0
                        pos = 0
                    # Initialise gap dict with 0 gaps
                    if rec['chrom'] not in gapdict:
                        gapdict[rec['chrom']] = 0
                    # If there is a gap in the alignment, record the size of the gap in gapdict
                    if int(rec['pos']) > pos:
                        # Add the gap size to gap dict
                        gapdict[rec['chrom']] += rec['pos'] - pos
                        # Set the expected position to the current position
                        pos = int(rec['pos'])
                    # Increment pos in preparation for the next iteration
                    pos += 1
                    # Initialise snpdict if necessary
                    if rec['chrom'] not in snpdict:
                        snpdict[rec['chrom']] = 0
                    # Initialise the current gene/allele in depthdict with the depth (reads_all) if necessary,
                    # otherwise add the current depth to the running total
                    if rec['chrom'] not in depthdict:
                        depthdict[rec['chrom']] = int(rec['reads_all'])
                    else:
                        depthdict[rec['chrom']] += int(rec['reads_all'])
                    # Dictionary of bases and the number of times each base was observed per position
                    bases = {'A': rec['A'], 'C': rec['C'], 'G': rec['G'], 'T': rec['T']}
                    # If the most prevalent base (calculated with max() and operator.itemgetter()) does not match the
                    # reference base, add this prevalent base to seqdict
                    if max(bases.items(), key=operator.itemgetter(1))[0] != rec['ref']:
                        seqdict[rec['chrom']] += max(bases.items(), key=operator.itemgetter(1))[0]
                        # Increment the running total of the number of SNPs
                        snpdict[rec['chrom']] += 1
                    else:
                        # If the bases match, add the reference base to seqdict
                        seqdict[rec['chrom']] += (rec['ref'])
                        # Initialise posdict if necessary, otherwise, increment the running total of matches
                        if rec['chrom'] not in matchdict:
                            matchdict[rec['chrom']] = 1
                        else:
                            matchdict[rec['chrom']] += 1
                    # Find the max and min coverage for each strain/gene combo
                    try:
                        maxdict[rec['chrom']] = int(rec['reads_all']) if \
                            int(rec['reads_all']) >= maxdict[rec['chrom']] else maxdict[rec['chrom']]
                    except KeyError:
                        maxdict[rec['chrom']] = int(rec['reads_all'])
                    try:
                        mindict[rec['chrom']] = int(rec['reads_all']) if \
                            int(rec['reads_all']) <= mindict[rec['chrom']] else mindict[rec['chrom']]
                    except KeyError:
                        mindict[rec['chrom']] = int(rec['reads_all'])
                    # Create a list of all the depths in order to calculate the standard deviation
                    try:
                        deviationdict[rec['chrom']].append(int(rec['reads_all']))
                    except KeyError:
                        deviationdict[rec['chrom']] = list()
                        deviationdict[rec['chrom']].append(int(rec['reads_all']))
            # If there are no results in the bam file, then pass over the strain
            except ValueError:
                pass
            # Iterate through all the genes/alleles with results above
            for allele in sorted(matchdict):
                # If the length of the match is greater or equal to the length of the gene/allele (multiplied by the
                # cutoff value) as determined using faidx indexing, then proceed
                if matchdict[allele] >= sample[analysistype].faidict[allele] * self.cutoff:
                    # Calculate the average depth by dividing the total number of reads observed by the
                    # length of the gene
                    averagedepth = float(depthdict[allele]) / float(matchdict[allele])
                    percentidentity = float(matchdict[allele]) / float(sample[analysistype].faidict[allele]) * 100
                    # Only report a positive result if this average depth is greater than 10X
                    if averagedepth > 10:
                        # Populate resultsdict with the gene/allele name, the percent identity, and the average depth
                        sample[analysistype].results.update({allele: '{:.2f}'.format(percentidentity)})
                        sample[analysistype].avgdepth.update({allele: '{:.2f}'.format(averagedepth)})
                        # Add the SNP and gap results to dictionaries
                        sample[analysistype].resultssnp.update({allele: snpdict[allele]})
                        sample[analysistype].resultsgap.update({allele: gapdict[allele]})
                        sample[analysistype].sequences.update({allele: seqdict[allele]})
                        sample[analysistype].maxcoverage.update({allele: maxdict[allele]})
                        sample[analysistype].mincoverage.update({allele: mindict[allele]})
                        sample[analysistype]\
                            .standarddev.update({allele: '{:.2f}'.format(numpy.std(deviationdict[allele], ddof=1))})
            self.parsequeue.task_done()

    def postmapping(self, analysistype, metadata):
        """

        :param analysistype:
        :param metadata:
        """
        import operator

        for sample in metadata:
            if sample[analysistype].results:
                if len(sample[analysistype].results) == 1:
                    for classification, percentidentity in sample[analysistype].results.items():
                        sample[self.analysistype].phylogeny.append(classification.split(analysistype.split('_')[0])[0]
                                                                   .rstrip('_'))
                        print('good', sample.name, classification, sample[self.analysistype].phylogeny, analysistype)
                else:
                    multidict = {key: float(value) for key, value in sample[analysistype].results.items()}
                    best = sorted(multidict.items(), key=operator.itemgetter(1), reverse=True)
                    print('multi', sample.name, multidict, best)
                    print(analysistype, analysistype.split('_')[0], best[0][0].split(analysistype.split('_')[0])[0]
                          .rstrip('_'))
                    sample[self.analysistype].phylogeny.append(best[0][0].split(analysistype.split('_')[0])[0]
                                                               .rstrip('_'))
                    # print(best[0][0].split('_')[0])
            else:
                sample[self.analysistype].complete = True
                print('no results', sample.name, sample[self.analysistype].phylogeny)

    def reporting(self):
        printtime('Creating reports', self.start)
        for sample in self.runmetadata:
            print(sample.name, sample[self.analysistype].phylogeny)


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
            # Sippr(self, self.cutoff)
            ProbeSippr(self, self.cutoff)
            #
            # self.attributer()
            # Create the reports
            # self.reporter()
            # Print the metadata
            printer = MetadataPrinter(self)
            printer.printmetadata()
            quit()

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

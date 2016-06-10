#!/usr/bin/env python
from glob import glob
from subprocess import call
from threading import Thread

from Bio.Sequencing.Applications import *

from SPAdesPipeline.OLCspades.accessoryFunctions import *
from SPAdesPipeline.OLCspades.bowtie import *

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
__author__ = 'adamkoziol'


class Custom(object):
    def targets(self):
        # There is a relatively strict databasing scheme necessary for the custom targets. Eventually, there will
        # be a helper script to combine individual files into a properly formatted combined file
        try:
            self.baitfile = glob('{}*.fasta'.format(self.targetpath))[0]
        # If the fasta file is missing, raise a custom error
        except IndexError as e:
            # noinspection PyPropertyAccess
            e.args = ['Cannot find the combined fasta file in {}. Please note that the file must have a '
                      '.fasta extension'.format(self.targetpath)]
            raise
        # Create the hash file of the baitfile
        targetbase = self.baitfile.split('.')[0]
        self.hashcall = 'cd {} && mirabait -b {} -k 19 -K {}.mhs.gz'.format(self.targetpath, self.baitfile, targetbase)
        if not os.path.isfile(self.baitfile.split('.')[0] + '.mhs.gz'):
            call(self.hashcall, shell=True, stdout=self.devnull, stderr=self.devnull)
        # Ensure that the hash file was successfully created
        self.hashfile = targetbase + '.mhs.gz'
        assert os.path.isfile(self.hashfile), u'Hashfile could not be created for the combined target file {0!r:s}' \
            .format(self.baitfile)
        # Bait
        self.baiting()

    def baiting(self):
        # Perform baiting
        printtime('Performing kmer baiting of fastq files with {} targets'.format(self.analysistype), self.start)
        # Create and start threads for each fasta file in the list
        for i in range(len(self.metadata)):
            # Send the threads to the bait method
            threads = Thread(target=self.bait, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Add the sample to the queue
                self.baitqueue.put(sample)
        self.baitqueue.join()
        # Run the bowtie2 read mapping module
        self.mapping()

    def bait(self):
        """
        Runs mirabait on the fastq files
        """
        while True:
            sample = self.baitqueue.get()
            # Create the .custom attribute
            setattr(sample, self.analysistype, GenObject())
            # Set attribute values
            sample[self.analysistype].outputdir = sample.general.outputdirectory + '/' + self.analysistype
            # Create the folder (if necessary)
            make_path(sample[self.analysistype].outputdir)
            # Set more attributes
            sample[self.analysistype].baitfile = self.baitfile
            sample[self.analysistype].hashfile = self.hashfile
            sample[self.analysistype].hashcall = self.hashcall
            sample[self.analysistype].targetpath = self.targetpath
            sample[self.analysistype].baitedfastq = '{}/{}_targetMatches.fastq'.format(sample[self.analysistype]
                                                                                       .outputdir, self.analysistype)

            # Make the system call
            if len(sample.general.fastqfiles) == 2:
                sample[self.analysistype].mirabaitcall = 'mirabait -c -B {} -t 4 -o {} -p {} {}' \
                    .format(sample[self.analysistype].hashfile, sample[self.analysistype].baitedfastq,
                            sample.general.fastqfiles[0], sample.general.fastqfiles[1])
            else:
                sample[self.analysistype].mirabaitcall = 'mirabait -c -B {} -t 4 -o {} {}' \
                    .format(sample[self.analysistype].hashfile, sample[self.analysistype].baitedfastq,
                            sample.general.fastqfiles[0])
            # Run the system call (if necessary)
            if not os.path.isfile(sample[self.analysistype].baitedfastq):
                call(sample[self.analysistype].mirabaitcall, shell=True, stdout=self.devnull, stderr=self.devnull)
            self.baitqueue.task_done()

    def mapping(self):
        printtime('Performing reference mapping', self.start)
        for i in range(len(self.metadata)):
            # Send the threads to
            threads = Thread(target=self.map, args=())
            # Set the daemon to True - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Set the path/name for the sorted bam file to be created
                sample[self.analysistype].sortedbam = '{}/{}_sorted.bam'.format(sample[self.analysistype].outputdir,
                                                                                self.analysistype)
                # Remove the file extension of the bait file for use in the indexing command
                sample[self.analysistype].baitfilenoext = sample[self.analysistype].baitfile.split('.')[0]
                # Use bowtie2 wrapper to create index the target file
                bowtie2build = Bowtie2BuildCommandLine(reference=sample[self.analysistype].baitfile,
                                                       bt2=sample[self.analysistype].baitfilenoext,
                                                       **self.builddict)
                # Use samtools wrapper to set up the bam sorting command
                samsort = SamtoolsSortCommandline(input_bam=sample[self.analysistype].sortedbam,
                                                  o=True,
                                                  out_prefix="-")
                # Create a list of programs to which data are piped as part of the reference mapping
                samtools = [
                    # When bowtie2 maps reads to all possible locations rather than choosing a 'best' placement, the
                    # SAM header for that read is set to 'secondary alignment', or 256. Please see:
                    # http://davetang.org/muse/2014/03/06/understanding-bam-flags/ The script below reads in the stdin
                    # and subtracts 256 from headers which include 256
                    'python {}/editsamheaders.py'.format(self.homepath),
                    # # Use samtools wrapper to set up the samtools view
                    SamtoolsViewCommandline(b=True,
                                            S=True,
                                            h=True,
                                            input_file="-"),
                    samsort]
                # Add custom parameters to a dictionary to be used in the bowtie2 alignment wrapper
                indict = {'--very-sensitive-local': True,
                          # For short targets, the match bonus can be increased
                          '--ma': self.matchbonus,
                          '-U': sample[self.analysistype].baitedfastq,
                          '-a': True,
                          '--threads': self.cpus,
                          '--local': True}
                # Create the bowtie2 reference mapping command
                bowtie2align = Bowtie2CommandLine(bt2=sample[self.analysistype].baitfilenoext,
                                                  threads=self.cpus,
                                                  samtools=samtools,
                                                  **indict)
                # Create the command to faidx index the bait file
                sample[self.analysistype].faifile = sample[self.analysistype].baitfile + '.fai'
                # In methods with multiple .fai files e.g. pathotyping different genera, this will be treated
                # differently
                self.faifile = sample[self.analysistype].faifile
                samindex = SamtoolsFaidxCommandline(reference=sample[self.analysistype].baitfile)
                # Add the commands (as strings) to the metadata
                sample[self.analysistype].bowtie2align = str(bowtie2align)
                sample[self.analysistype].bowtie2build = str(bowtie2build)
                sample[self.analysistype].samindex = str(samindex)
                # Add the commands to the queue. Note that the commands would usually be set as attributes of the sample
                # but there was an issue with their serialization when printing out the metadata
                if not os.path.isfile(sample[self.analysistype].baitfilenoext + '.1' + self.bowtiebuildextension):
                    stdoutbowtieindex, stderrbowtieindex = map(StringIO,
                                                               bowtie2build(cwd=sample[self.analysistype].targetpath))
                    # Write any error to a log file
                    if stderrbowtieindex:
                        # Write the standard error to log, bowtie2 puts alignment summary here
                        with open(os.path.join(sample[self.analysistype].targetpath,
                                               '{}_bowtie_index.log'.format(self.analysistype)), 'ab+') as log:
                            log.writelines(logstr(bowtie2build, stderrbowtieindex.getvalue(),
                                                  stdoutbowtieindex.getvalue()))
                    # Close the stdout and stderr streams
                    stdoutbowtieindex.close()
                    stderrbowtieindex.close()
                self.mapqueue.put((sample, bowtie2build, bowtie2align, samindex))
        self.mapqueue.join()
        # Use samtools to index the sorted bam file
        self.indexing()

    def map(self):
        while True:
            # Get the necessary values from the queue
            sample, bowtie2build, bowtie2align, samindex = self.mapqueue.get()
            # Use samtools faidx to index the bait file - this will be used in the sample parsing
            if not os.path.isfile(sample[self.analysistype].faifile):
                stdoutindex, stderrindex = map(StringIO, samindex(cwd=sample[self.analysistype].targetpath))
                # Write any error to a log file
                if stderrindex:
                    # Write the standard error to log, bowtie2 puts alignment summary here
                    with open(os.path.join(sample[self.analysistype].targetpath,
                                           '{}_samtools_index.log'.format(self.analysistype)), 'ab+') as log:
                        log.writelines(logstr(samindex, stderrindex.getvalue(), stdoutindex.getvalue()))
                # Close the stdout and stderr streams
                stdoutindex.close()
                stderrindex.close()
            # Only run the functions if the sorted bam files and the indexed bait file do not exist
            if not os.path.isfile(sample[self.analysistype].sortedbam):
                # Set stdout to a stringIO stream
                stdout, stderr = map(StringIO, bowtie2align(cwd=sample[self.analysistype].outputdir))
                if stderr:
                    # Write the standard error to log, bowtie2 puts alignment summary here
                    with open(os.path.join(sample[self.analysistype].outputdir,
                                           '{}_bowtie_samtools.log'.format(self.analysistype)), 'ab+') as log:
                        log.writelines(logstr(bowtie2align, stderr.getvalue(), stdout.getvalue()))
                stdout.close()
                stderr.close()
            self.mapqueue.task_done()

    def indexing(self):
        printtime('Indexing sorted bam files', self.start)
        for i in range(len(self.metadata)):
            # Send the threads to
            threads = Thread(target=self.index, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                bamindex = SamtoolsIndexCommandline(input=sample[self.analysistype].sortedbam)
                sample[self.analysistype].sortedbai = sample[self.analysistype].sortedbam + '.bai'
                sample[self.analysistype].bamindex = str(bamindex)
                self.indexqueue.put((sample, bamindex))
        self.indexqueue.join()
        # Parse the results
        self.parsing()

    def index(self):
        while True:
            sample, bamindex = self.indexqueue.get()
            # Only make the call if the .bai file doesn't already exist
            if not os.path.isfile(sample[self.analysistype].sortedbai):
                # Use cStringIO streams to handle bowtie output
                stdout, stderr = map(StringIO, bamindex(cwd=sample[self.analysistype].outputdir))
                if stderr:
                    # Write the standard error to log
                    with open(os.path.join(sample[self.analysistype].outputdir,
                                           '{}_samtools_bam_index.log'.format(self.analysistype)), 'ab+') as log:
                        log.writelines(logstr(bamindex, stderr.getvalue(), stdout.getvalue()))
                stderr.close()
            self.indexqueue.task_done()

    def parsing(self):
        printtime('Parsing sorted bam files', self.start)
        for i in range(len(self.metadata)):
            # Send the threads to
            threads = Thread(target=self.parse, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        # Get the fai file into a dictionary to be used in parsing results
        with open(self.faifile, 'rb') as faifile:
            for line in faifile:
                data = line.split('\t')
                self.faidict[data[0]] = int(data[1])
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                self.parsequeue.put(sample)
        self.parsequeue.join()

    def parse(self):
        import pysamstats
        import operator
        while True:
            sample = self.parsequeue.get()
            # Initialise dictionaries to store parsed data
            matchdict = dict()
            depthdict = dict()
            seqdict = dict()
            resultsdict = dict()
            snpdict = dict()
            gapdict = dict()
            snpresults = dict()
            gapresults = dict()
            seqresults = dict()
            # Variable to store the expected position in gene/allele
            pos = 0
            try:
                # Use the stat_variation function of pysam stats to return records parsed from sorted bam files
                # Values of interest can be retrieved using the appropriate keys
                for rec in pysamstats.stat_variation(alignmentfile=sample[self.analysistype].sortedbam,
                                                     fafile=sample[self.analysistype].baitfile,
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
                    # In strains with either multiple copies of the same gene, or multiple alleles with high sequence
                    # identity, the reference mapper cannot determine where to map the reads, and will map both to each
                    # sequence. This can be detected by finding positions where there are two bases that are almost
                    # equally represented in the pileup.
                    # TODO Implement variant calling
                    # Make a shallow copy of bases
                    # secondbases = dict(bases)
                    # Remove the most represented base in this copy
                    # del secondbases[max(bases.iteritems(), key=operator.itemgetter(1))[0]]
                    # If there are at least 30% as many reads with the second most prevalent bases compared to the
                    # most prevalent base
                    # if float(max(secondbases.iteritems(), key=operator.itemgetter(1))[1]) / \
                    #         float(max(bases.iteritems(), key=operator.itemgetter(1))[1]) > 0.3:
                    # VARIANTS!
                    #     print rec['pos'], bases, max(bases.iteritems(), key=operator.itemgetter(1)), \
                    #         max(secondbases.iteritems(), key=operator.itemgetter(1)), rec['ref']
                    # If the most prevalent base (calculated with max() and operator.itemgetter()) does not match the
                    # reference base, add this prevalent base to seqdict
                    if max(bases.iteritems(), key=operator.itemgetter(1))[0] != rec['ref']:
                        seqdict[rec['chrom']] += max(bases.iteritems(), key=operator.itemgetter(1))[0]
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
            # If there are no results in the bam file, then pass over the strain
            except ValueError:
                pass
            # Iterate through all the genes/alleles with results above
            for allele in sorted(matchdict):
                # If the length of the match is greater or equal to the length of the gene/allele (multiplied by the
                # cutoff value) as determined using faidx indexing, then proceed
                if matchdict[allele] >= self.faidict[allele] * self.cutoff:
                    # Calculate the average depth by dividing the total number of reads observed by the
                    # length of the gene
                    averagedepth = float(depthdict[allele]) / float(matchdict[allele])
                    percentidentity = float(matchdict[allele]) / float(self.faidict[allele]) * 100
                    # Only report a positive result if this average depth is greater than 4X
                    if averagedepth > 4:
                        # Populate resultsdict with the gene/allele name, the percent identity, and the average depth
                        resultsdict.update({allele: {'{:.2f}'.format(percentidentity): '{:.2f}'.format(averagedepth)}})
                        # Add the SNP and gap results to dictionaries
                        snpresults.update({allele: snpdict[allele]})
                        gapresults.update({allele: gapdict[allele]})
                        seqresults.update({allele: seqdict[allele]})
            # Add these results to the sample object
            sample[self.analysistype].results = resultsdict
            sample[self.analysistype].resultssnp = snpresults
            sample[self.analysistype].resultsgap = gapresults
            sample[self.analysistype].sequences = seqresults
            self.parsequeue.task_done()

    # noinspection PyDefaultArgument
    def __init__(self, inputobject, analysistype, cutoff=0.98, matchbonus=2, builddict=dict(), extension='.bt2'):
        from Queue import Queue
        from threading import Lock
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        self.targetpath = inputobject.customtargetpath if inputobject.customtargetpath else inputobject.targetpath
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.analysistype = analysistype
        self.cpus = inputobject.cpus
        self.homepath = inputobject.homepath
        self.baitfile = ''
        self.hashfile = ''
        self.hashcall = ''
        self.faifile = ''
        self.faidict = dict()
        self.devnull = open(os.devnull, 'wb')  # define /dev/null
        self.baitqueue = Queue(maxsize=self.cpus)
        self.mapqueue = Queue(maxsize=self.cpus)
        self.indexqueue = Queue(maxsize=self.cpus)
        self.parsequeue = Queue(maxsize=self.cpus)
        self.threadlock = Lock()
        self.cutoff = cutoff
        self.matchbonus = matchbonus
        self.builddict = builddict
        self.bowtiebuildextension = extension
        printtime('Performing analysis with {} targets folder'.format(self.analysistype), self.start)
        self.targets()

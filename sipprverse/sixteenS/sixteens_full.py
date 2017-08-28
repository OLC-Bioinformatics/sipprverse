#!/usr/bin/env python 3
from Bio import SeqIO
import operator
import time
from sipprCommon.sippingmethods import *
from sipprCommon.objectprep import Objectprep
from accessoryFunctions.accessoryFunctions import *

__author__ = 'adamkoziol'


class SixteenSBait(Sippr):

    def targets(self):
        """
        Create the GenObject for the analysis type, create the hash file for baiting (if necessary)
        """
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
                    'Hashfile could not be created for the target file {0!r:s}'.format(
                        sample[self.analysistype].baitfile)
                sample[self.analysistype].outputdir = os.path.join(sample.run.outputdirectory, self.analysistype)
                sample[self.analysistype].baitedfastq = \
                    '{}/{}_targetMatches.fastq'.format(sample[self.analysistype].outputdir, self.analysistype)
                sample[self.analysistype].complete = False
        # Run the baiting method
        self.baiting()

    def baiting(self):
        """
        Perform baiting of FASTQ files with mirabait
        """
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


class SixteenSSipper(Sippr):

    def targets(self):
        """
        Using the data from the BLAST analyses, set the targets folder, and create the 'mapping file'. This is the
        genera-specific FASTA file that will be used for all the reference mapping; it replaces the 'bait file' in the
        code
        """
        printtime('Performing analysis with {} targets folder'.format(self.analysistype), self.start)
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                sample[self.analysistype].targetpath = \
                    os.path.join(self.targetpath, 'genera', sample[self.analysistype].genus, '')
                # There is a relatively strict databasing scheme necessary for the custom targets. Eventually,
                # there will be a helper script to combine individual files into a properly formatted combined file
                try:
                    sample[self.analysistype].mappingfile = glob('{}*.fa'
                                                                 .format(sample[self.analysistype].targetpath))[0]
                # If the fasta file is missing, raise a custom error
                except IndexError as e:
                    # noinspection PyPropertyAccess
                    e.args = ['Cannot find the combined fasta file in {}. Please note that the file must have a '
                              '.fasta extension'.format(sample[self.analysistype].targetpath)]
                    if os.path.isdir(sample[self.analysistype].targetpath):
                        raise
                    else:
                        sample.general.bestassemblyfile = 'NA'
        # Run the reference mapping using the mapping file
        self.mapping()

    def mapping(self):
        """
        Perform reference mapping - this method overrides the default method in sipping methods. The major difference
        is the use of 'mapping file' instead of 'bait file'
        """
        printtime('Performing reference mapping', self.start)
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.map, args=())
            # Set the daemon to True - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                # Set the path/name for the sorted bam file to be created
                sample[self.analysistype].sortedbam = '{}/{}_sorted.bam'.format(sample[self.analysistype].outputdir,
                                                                                self.analysistype)
                # Remove the file extension of the bait file for use in the indexing command
                sample[self.analysistype].mappingfilenoext = sample[self.analysistype].mappingfile.split('.')[0]
                # Use bowtie2 wrapper to create index the target file
                bowtie2build = Bowtie2BuildCommandLine(reference=sample[self.analysistype].mappingfile,
                                                       bt2=sample[self.analysistype].mappingfilenoext,
                                                       **self.builddict)
                # Use samtools wrapper to set up the bam sorting command
                samsort = SamtoolsSortCommandline(input=sample[self.analysistype].sortedbam,
                                                  o=True,
                                                  out_prefix="-")
                # Determine the location of the SAM header editing script
                import sipprCommon.editsamheaders
                scriptlocation = sipprCommon.editsamheaders.__file__
                samtools = [
                    # When bowtie2 maps reads to all possible locations rather than choosing a 'best' placement, the
                    # SAM header for that read is set to 'secondary alignment', or 256. Please see:
                    # http://davetang.org/muse/2014/03/06/understanding-bam-flags/ The script below reads in the stdin
                    # and subtracts 256 from headers which include 256
                    'python3 {}'.format(scriptlocation),
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
                          '-U': sample[self.analysistype].baitedfastq,
                          '-a': True,
                          '--threads': self.threads,
                          '--local': True}
                # Create the bowtie2 reference mapping command
                bowtie2align = Bowtie2CommandLine(bt2=sample[self.analysistype].mappingfilenoext,
                                                  threads=self.threads,
                                                  samtools=samtools,
                                                  **indict)
                # Create the command to faidx index the bait file
                sample[self.analysistype].faifile = sample[self.analysistype].mappingfile + '.fai'
                samindex = SamtoolsFaidxCommandline(reference=sample[self.analysistype].mappingfile)
                # Add the samindex command (as a string) to the metadata - other commands don't seem to work properly,
                # and have been omitted
                sample[self.analysistype].samindex = str(samindex)
                # Add the commands to the queue. Note that the commands would usually be set as attributes of the sample
                # but there was an issue with their serialization when printing out the metadata
                if not os.path.isfile(sample[self.analysistype].mappingfilenoext + '.1' + self.bowtiebuildextension):
                    stdoutbowtieindex, stderrbowtieindex = map(StringIO,
                                                               bowtie2build(cwd=sample[self.analysistype].targetpath))
                    # Write any error to a log file
                    if stderrbowtieindex:
                        # Write the standard error to log, bowtie2 puts alignment summary here
                        with open(os.path.join(sample[self.analysistype].targetpath,
                                               '{}_bowtie_index.log'.format(self.analysistype)), 'a+') as log:
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
                                           '{}_samtools_index.log'.format(self.analysistype)), 'a+') as log:
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
                                           '{}_bowtie_samtools.log'.format(self.analysistype)), 'a+') as log:
                        log.writelines(logstr([bowtie2align], stderr.getvalue(), stdout.getvalue()))
                stdout.close()
                stderr.close()
            self.mapqueue.task_done()

    def parse(self):
        import pysamstats
        import numpy
        while True:
            sample = self.parsequeue.get()
            # Initialise dictionaries to store parsed data
            matchdict = dict()
            depthdict = dict()
            seqdict = dict()
            snpdict = dict()
            gapdict = dict()
            maxdict = dict()
            mindict = dict()
            deviationdict = dict()
            sample[self.analysistype].results = dict()
            sample[self.analysistype].avgdepth = dict()
            sample[self.analysistype].resultssnp = dict()
            sample[self.analysistype].resultsgap = dict()
            sample[self.analysistype].sequences = dict()
            sample[self.analysistype].maxcoverage = dict()
            sample[self.analysistype].mincoverage = dict()
            sample[self.analysistype].standarddev = dict()
            # Variable to store the expected position in gene/allele
            pos = 0
            try:
                # Use the stat_variation function of pysam stats to return records parsed from sorted bam files
                # Values of interest can be retrieved using the appropriate keys
                for rec in pysamstats.stat_variation(alignmentfile=sample[self.analysistype].sortedbam,
                                                     fafile=sample[self.analysistype].mappingfile,
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
                if matchdict[allele] >= sample[self.analysistype].faidict[allele] * self.cutoff:
                    # Calculate the average depth by dividing the total number of reads observed by the
                    # length of the gene
                    averagedepth = float(depthdict[allele]) / float(matchdict[allele])
                    percentidentity = float(matchdict[allele]) / float(sample[self.analysistype].faidict[allele]) * 100
                    # Only report a positive result if this average depth is greater than 10X
                    if averagedepth > 10:
                        # Populate resultsdict with the gene/allele name, the percent identity, and the average depth
                        sample[self.analysistype].results.update({allele: '{:.2f}'.format(percentidentity)})
                        sample[self.analysistype].avgdepth.update({allele: '{:.2f}'.format(averagedepth)})
                        # Add the SNP and gap results to dictionaries
                        sample[self.analysistype].resultssnp.update({allele: snpdict[allele]})
                        sample[self.analysistype].resultsgap.update({allele: gapdict[allele]})
                        sample[self.analysistype].sequences.update({allele: seqdict[allele]})
                        sample[self.analysistype].maxcoverage.update({allele: maxdict[allele]})
                        sample[self.analysistype].mincoverage.update({allele: mindict[allele]})
                        sample[self.analysistype]\
                            .standarddev.update({allele: '{:.2f}'.format(numpy.std(deviationdict[allele], ddof=1))})
            self.parsequeue.task_done()


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
            self.threads = int(self.cpus / len(self.runmetadata.samples)) \
                if self.cpus / len(self.runmetadata.samples) > 1 \
                else 1
            # Use a custom sippr method to use the full reference database as bait, and run mirabait against the FASTQ
            # reads - do not perform reference mapping yet
            SixteenSBait(self, self.cutoff)
            # Subsample 1000 reads from the FASTQ files
            self.subsample()
            # Convert the subsampled FASTQ files to FASTA format
            self.fasta()
            # Create BLAST databases if required
            self.makeblastdb()
            # Run BLAST analyses of the subsampled FASTA files against the NCBI 16S reference database
            self.blast()
            # Parse the BLAST results
            self.blastparse()
            # Feed the BLAST results into a modified sippr method to perform reference mapping using the calculated
            # genus of the sample as the mapping file
            SixteenSSipper(self, self.cutoff)
            # Create reports
            self.reporter()

    def subsample(self):
        """
        Subsample 1000 reads from the baited files
        """
        # Create the threads for the analysis
        printtime('Subsampling FASTQ reads', self.starttime)
        for _ in range(self.cpus):
            threads = Thread(target=self.subsamplethreads, args=())
            threads.setDaemon(True)
            threads.start()
        for sample in self.runmetadata.samples:
            # Set the name of the subsampled FASTQ file
            sample[self.analysistype].subsampledfastq = os.path.splitext(sample[self.analysistype].baitedfastq)[0] \
                                                        + '_subsampled.fastq'
            # Set the system call
            sample[self.analysistype].seqtkcall = 'seqtk sample {} 1000 > {}'\
                .format(sample[self.analysistype].baitedfastq,
                        sample[self.analysistype].subsampledfastq)
            # Add the sample to the queue
            self.samplequeue.put(sample)
        self.samplequeue.join()

    def subsamplethreads(self):
        while True:
            sample = self.samplequeue.get()
            # Check to see if the subsampled FASTQ file has already been created
            if not os.path.isfile(sample[self.analysistype].subsampledfastq):
                # Run the system call
                call(sample[self.analysistype].seqtkcall, shell=True, stdout=self.devnull, stderr=self.devnull)
            self.samplequeue.task_done()

    def fasta(self):
        """
        Convert the subsampled reads to FASTA format
        """
        printtime('Converting FASTQ files to FASTA format', self.starttime)
        # Create the threads for the analysis
        for _ in range(self.cpus):
            threads = Thread(target=self.fastathreads, args=())
            threads.setDaemon(True)
            threads.start()
        for sample in self.runmetadata.samples:
            # Set the name as the FASTA file - the same as the FASTQ, but with .fa file extension instead of .fastq
            sample[self.analysistype].fasta = os.path.splitext(sample[self.analysistype].subsampledfastq)[0] + '.fa'
            # Set the system call
            sample[self.analysistype].fastxcall = 'fastq_to_fasta -i {} -o {}'\
                .format(sample[self.analysistype].subsampledfastq, sample[self.analysistype].fasta)
            # Add the sample to the queue
            self.fastaqueue.put(sample)
        self.fastaqueue.join()

    def fastathreads(self):
        while True:
            sample = self.fastaqueue.get()
            # Check to see if the FASTA file already exists
            if not os.path.isfile(sample[self.analysistype].fasta):
                # Run the system call
                call(sample[self.analysistype].fastxcall, shell=True, stdout=self.devnull, stderr=self.devnull)
            self.fastaqueue.task_done()

    def makeblastdb(self):
        """
        Makes blast database files from targets as necessary
        """
        # Iterate through the samples to set the bait file. Break after the first entry, as all entries will be the same
        for sample in self.runmetadata.samples:
            self.baitfile = sample[self.analysistype].baitfile
            break
        # Remove the file extension
        db = self.baitfile.split('.')[0]
        # Add '.nhr' for searching below
        nhr = '{}.nhr'.format(db)
        # Check for already existing database files
        if not os.path.isfile(str(nhr)):
            # Create the databases
            call('makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'.format(self.baitfile, db),
                 shell=True, stdout=self.devnull, stderr=self.devnull)

    def blast(self):
        """
        Run BLAST analyses of the subsampled FASTQ reads against the NCBI 16S reference database
        """
        from Bio.Blast.Applications import NcbiblastnCommandline
        printtime('BLASTing FASTA files against {} database'.format(self.analysistype), self.starttime)
        for _ in range(self.cpus):
            threads = Thread(target=self.blastthreads, args=())
            threads.setDaemon(True)
            threads.start()
        for sample in self.runmetadata.samples:
            # Set the name of the BLAST report
            sample[self.analysistype].blastreport = os.path.join(
                sample[self.analysistype].outputdir, '{}_{}_blastresults.csv'.format(sample.name, self.analysistype))
            # Use the NCBI BLASTn command line wrapper module from BioPython to set the parameters of the search
            blastn = NcbiblastnCommandline(query=sample[self.analysistype].fasta,
                                           db=os.path.splitext(sample[self.analysistype].baitfile)[0],
                                           max_target_seqs=1,
                                           num_threads=12,
                                           outfmt="'6 qseqid sseqid positive mismatch gaps "
                                                  "evalue bitscore slen length qstart qend qseq sstart send sseq'",
                                           out=sample[self.analysistype].blastreport)
            # Add a string of the command to the metadata object
            sample[self.analysistype].blastcall = str(blastn)
            # Add the object and the command to the BLAST queue
            self.blastqueue.put((sample, blastn))
        self.blastqueue.join()

    def blastthreads(self):
        while True:
            sample, blastn = self.blastqueue.get()
            if not os.path.isfile(sample[self.analysistype].blastreport):
                # Ensure that the query file exists; this can happen with very small .fastq files
                if os.path.isfile(sample[self.analysistype].fasta):
                    # Perform the BLAST analysis
                    blastn()
            self.blastqueue.task_done()

    def blastparse(self):
        """
        Parse the blast results, and store necessary data in dictionaries in sample object
        """
        printtime('Parsing BLAST results', self.starttime)
        from csv import DictReader
        # Load the NCBI 16S reference database as a dictionary
        dbrecords = SeqIO.to_dict(SeqIO.parse(self.baitfile, 'fasta'))
        for sample in self.runmetadata.samples:
            # Allow for no BLAST results
            if os.path.isfile(sample[self.analysistype].blastreport):
                # Initialise a dictionary to store the number of times a genus is the best hit
                sample[self.analysistype].frequency = dict()
                # Open the sequence profile file as a dictionary
                blastdict = DictReader(open(sample[self.analysistype].blastreport),
                                       fieldnames=self.fieldnames, dialect='excel-tab')
                for record in blastdict:
                    # Create the subject id. It will look like this: gi|1018196593|ref|NR_136472.1|
                    subject = record['subject_id']
                    # Extract the genus name. Use the subject id as a key in the dictionary of the reference database.
                    # It will return the full record e.g. gi|1018196593|ref|NR_136472.1| Escherichia marmotae
                    # strain HT073016 16S ribosomal RNA, partial sequence
                    # This full description can be manipulated to extract the genus e.g. Escherichia
                    genus = dbrecords[subject].description.split('|')[-1].split()[0]
                    # Increment the number of times this genus was encountered, or initialise the dictionary with this
                    # genus the first time it is seen
                    try:
                        sample[self.analysistype].frequency[genus] += 1
                    except KeyError:
                        sample[self.analysistype].frequency[genus] = 1
                # Sort the dictionary based on the number of times a genus is seen
                sample[self.analysistype].sortedgenera = sorted(sample[self.analysistype].frequency.items(),
                                                                key=operator.itemgetter(1), reverse=True)
                # Extract the top result, and set it as the genus of the sample
                sample[self.analysistype].genus = sample[self.analysistype].sortedgenera[0][0]
                # Previous code relies on having the closest refseq genus, so set this as above
                sample.general.closestrefseqgenus = sample[self.analysistype].genus
            else:
                # Populate attributes with 'NA'
                sample[self.analysistype].sortedgenera = 'NA'
                sample[self.analysistype].genus = 'NA'
                sample.general.closestrefseqgenus = 'NA'

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # Initialise the header and data strings
        header = 'Strain,Gene,PercentIdentity,Genus,FoldCoverage\n'
        data = ''
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                try:
                    # Select the best hit of all the full-length 16S genes mapped
                    sample[self.analysistype].besthit = sorted(sample[self.analysistype].results.items(),
                                                               key=operator.itemgetter(1), reverse=True)[0][0]
                    # Add the sample name to the data string
                    data += sample.name + ','
                    # Find the record that matches the best hit, and extract the necessary values to be place in the
                    # data string
                    for name, identity in sample[self.analysistype].results.items():
                        if name == sample[self.analysistype].besthit:
                            data += '{},{},{},{}\n'.format(name, identity, sample[self.analysistype].genus,
                                                           sample[self.analysistype].avgdepth[name])
                except (KeyError, IndexError):
                    data += '\n'
            # Write the results to the report
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
        from queue import Queue
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
        self.threads = int()
        self.runmetadata = args.runmetadata
        self.pipeline = args.pipeline
        self.copy = args.copy
        self.devnull = open(os.path.devnull, 'w')
        self.samplequeue = Queue(maxsize=self.cpus)
        self.fastaqueue = Queue(maxsize=self.cpus)
        self.blastqueue = Queue(maxsize=self.cpus)
        self.baitfile = str()
        self.taxonomy = {'Escherichia': 'coli', 'Listeria': 'monocytogenes', 'Salmonella': 'enterica'}
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'query_sequence',
                           'subject_start', 'subject_end', 'subject_sequence']
        # Run the analyses
        self.runner()

if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    from subprocess import Popen
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                   shell=True, stdout=PIPE).communicate()[0].rstrip()
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

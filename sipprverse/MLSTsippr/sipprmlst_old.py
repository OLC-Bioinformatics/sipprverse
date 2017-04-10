#!/usr/bin/env python
from SPAdesPipeline.OLCspades.mMLST import *
from subprocess import call
# from customtargets import *

__author__ = 'adamkoziol'


class MLSTmap(object):
    def targets(self):
        printtime('Finding {} target files'.format(self.analysistype), self.start)
        for sample in self.runmetadata:
            setattr(sample, self.analysistype, GenObject())
            if self.analysistype.lower() == 'rmlst':
                # Run the allele updater method
                updatecall, allelefolder = getrmlsthelper(self.referencefilepath, self.updatermlst, self.start)
                self.alleles = glob('{}/*.tfa'.format(allelefolder))
                self.profile = glob('{}/*.txt'.format(allelefolder))
                self.supplementalprofile = '{}rMLST/OLC_rMLST_profiles.txt'.format(self.referencefilepath)
                self.combinedalleles = glob('{}/*.fasta'.format(allelefolder))[0]
                # Set the metadata file appropriately
                sample[self.analysistype].alleledir = allelefolder
                sample[self.analysistype].updatecall = updatecall
                sample[self.analysistype].targetpath = allelefolder
            else:
                try:
                    allelefolder = sorted(glob('{}MLST/{}/*/'.format(self.referencefilepath,
                                                                     sample.mash.closestrefseqgenus)))[-1]
                except IndexError:
                    allelefolder = 'NA'
                self.alleles = glob('{}*.tfa'.format(allelefolder))
                self.profile = glob('{}*.txt'.format(allelefolder))
                try:
                    self.combinedalleles = glob('{}*.fasta'.format(allelefolder))[0]
                except IndexError:
                    # Set the metadata file appropriately
                    sample.general.bestassemblyfile = 'NA'
                sample[self.analysistype].alleledir = allelefolder
            # Add the combined alleles to the profile set
            self.profileset.add(self.combinedalleles)
            sample[self.analysistype].alleles = self.alleles
            sample[self.analysistype].allelenames = [os.path.split(x)[1].split('.')[0] for x in self.alleles]
            sample[self.analysistype].profile = self.profile
            sample[self.analysistype].analysistype = self.analysistype
            sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory, self.analysistype)
            sample[self.analysistype].combinedalleles = self.combinedalleles
            sample[self.analysistype].supplementalprofile = self.supplementalprofile
            sample[self.analysistype].baitfile = sample[self.analysistype].combinedalleles

        # Process the targets
        printtime('Indexing {} target file'.format(self.analysistype), self.start)
        for target in self.profileset:
            # Create the hash file of the baitfile
            targetbase = target.split('.')[0]
            hashcall = 'cd {} && mirabait -b {} -k 31 -K {}.mhs.gz'.format(self.targetpath, target, targetbase)
            hashfile = targetbase + '.mhs.gz'
            if not os.path.isfile(hashfile):
                call(hashcall, shell=True, stdout=self.devnull, stderr=self.devnull)
            # Ensure that the hash file was successfully created
            # assert os.path.isfile(hashfile), u'Hashfile could not be created for the combined target file {0!r:s}' \
            #     .format(target)
            # Populate the appropriate attributes
            for sample in self.runmetadata:
                if sample.general.bestassemblyfile != 'NA':
                    if sample[self.analysistype].combinedalleles == target:
                        sample[self.analysistype].hashcall = hashcall
                        sample[self.analysistype].hashfile = hashfile
        # Bait
        self.baiting()

    def baiting(self):
        # Perform baiting
        printtime('Performing kmer baiting of fastq files with targets', self.start)
        # Create and start threads for each fasta file in the list
        # for i in range(len(self.runmetadata)):
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
        # Run the bowtie2 read mapping module
        self.mapping()

    def bait(self):
        """
        Runs mirabait on the fastq files
        """
        while True:
            sample = self.baitqueue.get()
            # Set attribute values
            sample[self.analysistype].targetpath = self.targetpath
            sample[self.analysistype].outputdir = sample.general.outputdirectory + '/' + self.analysistype
            sample[self.analysistype].baitedfastq = '{}/{}_targetMatches.fastq'.format(sample[self.analysistype]
                                                                                       .outputdir, self.analysistype)
            # Create the folder (if necessary)
            make_path(sample[self.analysistype].outputdir)
            # Create the system call
            if len(sample.general.fastqfiles) == 2:
                sample[self.analysistype].mirabaitcall = 'mirabait -c -B {} -t 4 -m 2048 -o {} -p {} {}' \
                    .format(sample[self.analysistype].hashfile, sample[self.analysistype].baitedfastq,
                            sample.general.fastqfiles[0], sample.general.fastqfiles[1])
            else:
                sample[self.analysistype].mirabaitcall = 'mirabait -c -B {} -t 4 -m 2048 -o {} {}' \
                    .format(sample[self.analysistype].hashfile, sample[self.analysistype].baitedfastq,
                            sample.general.fastqfiles[0])
            # Run the system call (if necessary)
            if not os.path.isfile(sample[self.analysistype].baitedfastq):
                call(sample[self.analysistype].mirabaitcall, shell=True, stdout=self.devnull, stderr=self.devnull)
            self.baitqueue.task_done()

    def mapping(self):
        """
        Completely changed the mapping logic, so this method overrides the default method
        """
        self.rmlstsketching()
        self.rmlstmashing()

    def rmlstsketching(self):
        """
        Create a sketch file of the baited fastq to be used by mash to determine the closest alleles for each gene
        """
        printtime('Indexing {} sorted bam files'.format(self.analysistype), self.start)
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.rmlstsketch, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                sample[self.analysistype].sketchcall = 'mash sketch -s 100000 -m 2 -p {} -o {}{}_sketch {}'.format(
                    self.cpus, sample[self.analysistype].reportdir, self.analysistype,
                    sample[self.analysistype].baitedfastq)
                sample[self.analysistype].sketchfile = '{}{}_sketch.msh'.format(sample[self.analysistype].reportdir,
                                                                                self.analysistype)
                self.sketchqueue.put(sample)
        self.sketchqueue.join()

    def rmlstsketch(self):
        while True:
            sample = self.sketchqueue.get()
            if not os.path.isfile(sample[self.analysistype].sketchfile):
                call(sample[self.analysistype].sketchcall, shell=True, stdout=self.devnull, stderr=self.devnull)
            self.sketchqueue.task_done()

    def rmlstmashing(self):
        """
        Run mash to determine the closest alleles for each gene in the analysis
        """
        printtime('Finding closest alleles for each {} gene target'.format(self.analysistype), self.start)
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                if 'alleles' in sample[self.analysistype].datastore:
                    for i in range(len(sample[self.analysistype].allelenames)):
                        # Send the threads to
                        threads = Thread(target=self.rmlstmash, args=())
                        # Set the daemon to true - something to do with thread management
                        threads.setDaemon(True)
                        # Start the threading
                        threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                # Set the name, and create the directory to store the mash tables
                sample[self.analysistype].mashtabledir = os.path.join(sample[self.analysistype].reportdir, 'mashtables')
                make_path(sample[self.analysistype].mashtabledir)
                outtablelist = list()
                if 'alleles' in sample[self.analysistype].datastore:
                    for allele in sample[self.analysistype].alleles:
                        # Create the name for the table file
                        outtable = '{}/{}.tab'.format(sample[self.analysistype].mashtabledir,
                                                      os.path.split(allele)[1].split('.')[0])
                        outtablelist.append(outtable)
                        # Not adding the mash command to the object, as there are 53 rMLST genes, and it will get messy
                        mashcommand = 'mash dist -i -p {} {} {} | sort -gk3 > {}'.format(
                            self.cpus, sample[self.analysistype].sketchfile, allele, outtable)
                        self.mashqueue.put((outtable, mashcommand))
                    sample[self.analysistype].outtables = sorted(outtablelist)
        self.mashqueue.join()
        self.rmlstmashparsing()

    def rmlstmash(self):
        while True:
            outtable, mashcommand = self.mashqueue.get()
            # Run the command if the table has not yet been created
            if not os.path.isfile(outtable):
                call(mashcommand, shell=True, stdout=self.devnull, stderr=self.devnull)
            self.mashqueue.task_done()

    def rmlstmashparsing(self):
        """
        Parses mash results to determine the five closest alleles for each gene in the analysis
        """
        printtime('Parsing closest {} allele matches'.format(self.analysistype), self.start)
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.rmlstmashparse, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                sample[self.analysistype].mashalleles = list()
                if 'alleles' in sample[self.analysistype].datastore:
                    for table in sample[self.analysistype].outtables:
                        self.mashparsequeue.put((sample, table))
        self.mashparsequeue.join()
        self.reduceddatabasecreating()

    def rmlstmashparse(self):
        while True:
            sample, table = self.mashparsequeue.get()
            # Open the mash results and extract the top five lines
            data = open(table, 'rb').readlines()[:5]
            for row in data:
                # Populate the attribute with the gene/allele name from the mash results
                sample[self.analysistype].mashalleles.append(row.split('\t')[1])
            self.mashparsequeue.task_done()

    def reduceddatabasecreating(self):
        """
        Uses results from mash to create a database of the five closest alleles for each gene
        """
        printtime('Reducing {} database'.format(self.analysistype), self.start)
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.reduceddatabasecreator, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                sample[self.analysistype].reduceddatabase = '{}/{}_reduceddatabase.fasta'.format(
                    sample[self.analysistype].outputdir, self.analysistype)
                if 'alleles' in sample[self.analysistype].datastore:
                    self.databasesqueue.put(sample)
        self.databasesqueue.join()
        self.databaseindexing()

    def reduceddatabasecreator(self):
        from Bio import SeqIO
        while True:
            sample = self.databasesqueue.get()
            # Only create the reduced database if it doesn't already exist
            if not os.path.isfile(sample[self.analysistype].reduceddatabase):
                rmlstdatabase = SeqIO.parse(sample[self.analysistype].baitfile, 'fasta')
                SeqIO.write((allele for allele in rmlstdatabase if allele.id in sample[self.analysistype].mashalleles),
                            sample[self.analysistype].reduceddatabase, 'fasta')
            self.databasesqueue.task_done()

    def databaseindexing(self):
        printtime('Performing {} reference mapping'.format(self.analysistype), self.start)
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
                sample[self.analysistype].databasenoext = sample[self.analysistype].reduceddatabase.split('.')[0]
                # Use bowtie2 wrapper to create index the target file
                bowtie2build = Bowtie2BuildCommandLine(reference=sample[self.analysistype].reduceddatabase,
                                                       bt2=sample[self.analysistype].databasenoext,
                                                       **self.builddict)
                # Use samtools wrapper to set up the bam sorting command
                samsort = SamtoolsSortCommandline(input_bam=sample[self.analysistype].sortedbam,
                                                  o=True,
                                                  out_prefix="-")
                # Create a list of programs to which data are piped as part of the reference mapping
                samtools = [
                    # When bowtie2 maps reads to all possible locations rather than just choosing a "best" placement,
                    # the SAM header for that read is set to 'secondary alignment', or 256. Please see:
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
                          '--local': True}
                # Create the bowtie2 reference mapping command
                bowtie2align = Bowtie2CommandLine(bt2=sample[self.analysistype].databasenoext,
                                                  threads=self.cpus,
                                                  samtools=samtools,
                                                  **indict)
                # Create the command to faidx index the bait file
                sample[self.analysistype].faifile = sample[self.analysistype].reduceddatabase + '.fai'
                samindex = SamtoolsFaidxCommandline(reference=sample[self.analysistype].reduceddatabase)
                # Add the commands (as strings) to the metadata
                sample[self.analysistype].bowtie2align = str(bowtie2align)
                sample[self.analysistype].bowtie2build = str(bowtie2build)
                sample[self.analysistype].samindex = str(samindex)
                # Add the commands to the queue. Note that the commands would usually be set as attributes of the sample
                # but there was an issue with their serialization when printing out the metadata
                if not os.path.isfile(sample[self.analysistype].databasenoext + '.1' + self.bowtiebuildextension):
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

    def parsing(self):
        printtime('Parsing {} sorted bam files'.format(self.analysistype), self.start)
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.parse, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        # Get the fai file into a dictionary to be used in parsing results
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                sample[self.analysistype].faidict = dict()
                with open(sample[self.analysistype].faifile, 'rb') as faifile:
                    for line in faifile:
                        data = line.split('\t')
                        sample[self.analysistype].faidict[data[0]] = int(data[1])
                self.parsequeue.put(sample)
        self.parsequeue.join()
        self.profiler()

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
            genespresent = set()
            closematches = dict()
            # Variable to store the expected position in gene/allele
            pos = 0
            try:
                # Use the stat_variation function of pysam stats to return records parsed from sorted bam files
                # Values of interest can be retrieved using the appropriate keys
                for rec in pysamstats.stat_variation(alignmentfile=sample[self.analysistype].sortedbam,
                                                     fafile=sample[self.analysistype].reduceddatabase,
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
                    # If there is an indel in the alignment, record it in gapdict
                    if float(rec['insertions']) / float(rec['reads_all']) > 0.5:
                        gapdict[rec['chrom']] += 1
                    if float(rec['deletions']) / float(rec['reads_all']) > 0.5:
                        gapdict[rec['chrom']] += 1
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
                try:
                    # Calculate the average depth by dividing the total number of reads observed by the
                    # length of the gene and percent identity by dividing the length of the match by the length of
                    # the reference allele sequence
                    averagedepth = float(depthdict[allele]) / float(matchdict[allele])
                    percentidentity = float(matchdict[allele]) / float(sample[self.analysistype].faidict[allele]) * 100
                except KeyError:
                    pass
                # Only report a positive result if this average depth is greater than 4X
                if averagedepth > 4:
                    # If the sequence has a 100% identity, and there are no indels, proceed
                    if matchdict[allele] >= sample[self.analysistype].faidict[allele] and gapdict[allele] == 0:
                        # Populate resultsdict with the gene/allele name, the percent identity, and the average depth
                        allelename = allele.split('_')[0] if '_' in allele else allele.split('-')[0]
                        resultsdict.update({allele: {'{:.2f}'.format(percentidentity): '{:.2f}'.format(averagedepth)}})
                        genespresent.add(allelename)
                    elif matchdict[allele] >= sample[self.analysistype].faidict[allele] * self.cutoff:
                        closematches.update({allele: {'{:.2f}'.format(percentidentity): '{:.2f}'.format(averagedepth)}})
            # Initialise a copy of close matches to remove any genes that have exact matches
            filteredclosematches = dict(closematches)
            for gene in sorted(sample[self.analysistype].allelenames):
                if gene in genespresent:
                    for allele in sorted(closematches):
                        if gene in allele:
                            try:
                                del filteredclosematches[allele]
                            except KeyError:
                                pass

            for gene in sample[self.analysistype].allelenames:
                # Initialise a variable to determine whether all the genes in the database have at least one hit
                # foundallele = False
                # If there are no perfect matches
                if gene not in genespresent:
                    try:
                        # Sort the closest matches by percent identity, and choose the top result
                        allele = max(filteredclosematches.iteritems(), key=operator.itemgetter(1))[0]
                        # If the gene name is within the allele string
                        if gene in allele:
                            # The gene has a closest match
                            # foundallele = True
                            # Calculate the percent identity and average depth as above
                            percentidentity = float(max(closematches.iteritems(),
                                                        key=operator.itemgetter(1))[1].items()[0][0])
                            averagedepth = float(max(closematches.iteritems(),
                                                     key=operator.itemgetter(1))[1].items()[0][1])
                            # Update the results
                            resultsdict.update({allele: {'{:.2f}'.format(percentidentity): '{:.2f}'
                                               .format(averagedepth)}})
                            # Add the SNP and gap results to dictionaries
                            snpresults.update({allele: snpdict[allele]})
                            gapresults.update({allele: gapdict[allele]})
                            # Store the sequence of the observed sequence
                            seqresults.update({gene: seqdict[allele]})
                    # Populate the results with 'negative' values
                    except ValueError:
                        resultsdict.update({gene: {'N': 0}})
            # Add these results to the sample object
            sample[self.analysistype].results = resultsdict
            sample[self.analysistype].newalleles = seqresults
            sample[self.analysistype].resultssnp = snpresults
            sample[self.analysistype].resultsgap = gapresults
            self.parsequeue.task_done()

    def profiler(self):
        """Creates a dictionary from the profile scheme(s)"""
        printtime('Loading {} sequence profiles'.format(self.analysistype), self.start)
        # Initialise variables
        profiledata = defaultdict(make_dict)
        profileset = set()
        supplementalset = ''
        genedict = {}
        # Find all the unique profiles to use with a set
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                if sample[self.analysistype].profile != 'NA':
                    profileset.add(sample[self.analysistype].profile[0])
                    if self.analysistype == 'rmlst':
                        supplementalset = sample[self.analysistype].supplementalprofile

        # Extract the profiles for each set
        for sequenceprofile in profileset:
            # Clear the list of genes
            genelist = []
            for sample in self.runmetadata:
                if sample.general.bestassemblyfile != 'NA':
                    if sequenceprofile == sample[self.analysistype].profile[0]:
                        genelist = [os.path.split(x)[1].split('.')[0] for x in sample[self.analysistype].alleles]
            try:
                # Open the sequence profile file as a dictionary
                profile = DictReader(open(sequenceprofile), dialect='excel-tab')
            # Revert to standard comma separated values
            except KeyError:
                # Open the sequence profile file as a dictionary
                profile = DictReader(open(sequenceprofile))
            # Iterate through the rows
            for row in profile:
                # Iterate through the genes
                for gene in genelist:
                    # Add the sequence profile, and type, the gene name and the allele number to the dictionary
                    try:
                        profiledata[sequenceprofile][row['ST']][gene] = row[gene]
                    except KeyError:
                        profiledata[sequenceprofile][row['rST']][gene] = row[gene]
            # Load the supplemental profile definitions
            if self.analysistype == 'rmlst':
                supplementalprofile = DictReader(open(supplementalset), dialect='excel-tab')
                # Do the same with the supplemental profile
                for row in supplementalprofile:
                    # Iterate through the genes
                    for gene in genelist:
                        # Add the sequence profile, and type, the gene name and the allele number to the dictionary
                        profiledata[sequenceprofile][row['rST']][gene] = row[gene]
            # Add the gene list to a dictionary
            genedict[sequenceprofile] = sorted(genelist)
            # Add the profile data, and gene list to each sample
            for sample in self.runmetadata:
                if sample.general.bestassemblyfile != 'NA':
                    if sequenceprofile == sample[self.analysistype].profile[0]:
                        # Populate the metadata with the profile data
                        sample[self.analysistype].profiledata = profiledata[sample[self.analysistype].profile[0]]
                        # Add the allele directory to a list of directories used in this analysis
                        self.allelefolders.add(sample[self.analysistype].alleledir)
                    dotter()
        self.sequencetyping()

    def sequencetyping(self):
        printtime('Determining {} sequence types'.format(self.analysistype), self.start)
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.sequencetyper, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                if 'results' in sample[self.analysistype].datastore:
                    self.typequeue.put(sample)
            # Populate the object with negative results
            else:
                sample[self.analysistype].sequencetype = 'NA'
                sample[self.analysistype].matchestosequencetype = 'NA'
                sample[self.analysistype].mismatchestosequencetype = 'NA'
        self.typequeue.join()
        # Run the report creation method
        self.reporter()

    def sequencetyper(self):
        while True:
            sample = self.typequeue.get()
            # Initialise variables
            header = 0
            # Iterate through the genomes
            genome = sample.name
            # Initialise self.bestmatch[genome] with an int that will eventually be replaced by the # of matches
            self.bestmatch[genome] = defaultdict(int)
            if sample[self.analysistype].profile != 'NA':
                # Create the profiledata variable to avoid writing self.profiledata[self.analysistype]
                # profiledata = self.profiledata[self.analysistype]
                profiledata = sample[self.analysistype].profiledata
                # For each gene name in the list of gene names
                for gene in sample[self.analysistype].allelenames:
                    # Clear the appropriate count and lists
                    multiallele = []
                    multipercent = []
                    # Go through the alleles
                    for geneallele in sample[self.analysistype].results:
                        if gene in geneallele:
                            try:
                                allele = geneallele.split('_')[1] if '_' in geneallele else geneallele.split('-')[1]
                            except IndexError:
                                allele = 'N'
                            percentid = sample[self.analysistype].results[geneallele].items()[0][0]
                            # "N" alleles screw up the allele splitter function
                            if allele != "N":
                                # Append as appropriate - alleleNumber is treated as an integer for proper sorting
                                multiallele.append(int(allele))
                                multipercent.append(percentid)
                            # If the allele is "N"
                            else:
                                # Append "N" and a percent identity of 0
                                multiallele.append("N")
                                multipercent.append(0)
                            if not multiallele:
                                multiallele.append("N")
                                multipercent.append(0)

                    # For whatever reason, the rMLST profile scheme treat multiple allele hits as 'N's.
                    # if len(multiallele) > 1:
                    #     print gene, sorted(multiallele)
                    multiallele = multiallele if len(multiallele) >= 1 else ['N']
                    if multipercent:
                        multipercent = multipercent if len(multiallele) == 1 else [0, 0]
                    else:
                        multipercent = [0]
                    # Populate self.bestdict with genome, gene, alleles joined with a space (this was made like
                    # this because allele is a list generated by the .iteritems() above
                    self.bestdict[genome][gene][" ".join(str(allele)
                                                         for allele in sorted(multiallele))] = multipercent[0]
                    # Find the profile with the most alleles in common with the query genome
                    for sequencetype in profiledata:
                        # The number of genes in the analysis
                        header = len(profiledata[sequencetype])
                        # refallele is the allele number of the sequence type
                        refallele = profiledata[sequencetype][gene]
                        # If there are multiple allele matches for a gene in the reference profile e.g. 10 692
                        if len(refallele.split(" ")) > 1:
                            # Map the split (on a space) alleles as integers - if they are treated as integers,
                            # the alleles will sort properly
                            intrefallele = map(int, refallele.split(" "))
                            # Create a string of the joined, sorted alleles
                            sortedrefallele = " ".join(str(allele) for allele in sorted(intrefallele))
                        else:
                            # Use the reference allele as the sortedRefAllele
                            sortedrefallele = refallele
                        for allele, percentid in self.bestdict[genome][gene].iteritems():
                            # If the allele in the query genome matches the allele in the reference profile, add
                            # the result to the bestmatch dictionary. Genes with multiple alleles were sorted
                            # the same, strings with multiple alleles will match: 10 692 will never be 692 10
                            if allele == sortedrefallele:
                                # Increment the number of matches to each profile
                                self.bestmatch[genome][sequencetype] += 1
                            # Special handling of BACT000060 and BACT000065 genes. When the reference profile
                            # has an allele of 'N', and the query allele doesn't, set the allele to 'N', and
                            # count it as a match
                            elif gene == 'BACT000060' or gene == 'BACT000065':
                                if sortedrefallele == 'N' and allele != 'N':
                                    # Increment the number of matches to each profile
                                    self.bestmatch[genome][sequencetype] += 1
                # Get the best number of matches
                # From: https://stackoverflow.com/questions/613183/sort-a-python-dictionary-by-value
                try:
                    sortedmatches = sorted(self.bestmatch[genome].items(), key=operator.itemgetter(1),
                                           reverse=True)[0][1]
                # If there are no matches, set :sortedmatches to zero
                except IndexError:
                    sortedmatches = 0
                # Otherwise, the query profile matches the reference profile
                if int(sortedmatches) == header:
                    # Iterate through best match
                    for sequencetype, matches in self.bestmatch[genome].iteritems():
                        if matches == sortedmatches:
                            for gene in profiledata[sequencetype]:
                                # Populate resultProfile with the genome, best match to profile, # of matches
                                # to the profile, gene, query allele(s), reference allele(s), and % identity
                                self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                    self.bestdict[genome][gene]
                                        .keys()[0]] = str(self.bestdict[genome][gene].values()[0])
                            sample[self.analysistype].sequencetype = sequencetype
                            sample[self.analysistype].matchestosequencetype = matches
                # If there are fewer matches than the total number of genes in the typing scheme
                elif 0 < int(sortedmatches) < header:
                    mismatches = []
                    # Iterate through the sequence types and the number of matches in bestDict for each genome
                    for sequencetype, matches in self.bestmatch[genome].iteritems():
                        # If the number of matches for a profile matches the best number of matches
                        if matches == sortedmatches:
                            # Iterate through the gene in the analysis
                            for gene in profiledata[sequencetype]:
                                # Get the reference allele as above
                                refallele = profiledata[sequencetype][gene]
                                # As above get the reference allele split and ordered as necessary
                                if len(refallele.split(" ")) > 1:
                                    intrefallele = map(int, refallele.split(" "))
                                    sortedrefallele = " ".join(str(allele) for allele in sorted(intrefallele))
                                else:
                                    sortedrefallele = refallele
                                # Populate self.mlstseqtype with the genome, best match to profile, # of matches
                                # to the profile, gene, query allele(s), reference allele(s), and % identity
                                self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                    self.bestdict[genome][gene].keys()[0]] \
                                    = str(self.bestdict[genome][gene].values()[0])
                                if sortedrefallele != self.bestdict[sample.name][gene].keys()[0]:
                                    mismatches.append(
                                        ({gene: ('{} ({})'.format(self.bestdict[sample.name][gene]
                                                                  .keys()[0], sortedrefallele))}))
                                sample[self.analysistype].mismatchestosequencetype = mismatches
                                sample[self.analysistype].sequencetype = sequencetype
                                sample[self.analysistype].matchestosequencetype = matches
                elif sortedmatches == 0:
                    for gene in sample[self.analysistype].allelenames:
                        # Populate the results profile with negative values for sequence type and sorted matches
                        self.resultprofile[genome]['NA'][sortedmatches][gene]['NA'] = 0
                    # Add the new profile to the profile file (if the option is enabled)
                    sample[self.analysistype].sequencetype = 'NA'
                    sample[self.analysistype].matchestosequencetype = 'NA'
                    sample[self.analysistype].mismatchestosequencetype = 'NA'
                dotter()
            self.typequeue.task_done()

    def reporter(self):
        """ Parse the results into a report"""
        printtime('Creating {} reports'.format(self.analysistype), self.start)
        # Initialise variables
        combinedrow = ''
        reportdirset = set()
        # Populate a set of all the report directories to use. A standard analysis will only have a single report
        # directory, while pipeline analyses will have as many report directories as there are assembled samples
        for sample in self.runmetadata:
            # Ignore samples that lack a populated reportdir attribute
            if sample[self.analysistype].reportdir != 'NA':
                make_path(sample[self.analysistype].reportdir)
                # Add to the set - I probably could have used a counter here, but I decided against it
                reportdirset.add(sample[self.analysistype].reportdir)
        # Create a report for each sample from :self.resultprofile
        for sample in self.runmetadata:
            if sample[self.analysistype].reportdir != 'NA':
                if type(sample[self.analysistype].allelenames) == list and sample.general.bestassemblyfile != 'NA':
                    # Populate the header with the appropriate data, including all the genes in the list of targets
                    row = 'Strain,Genus,SequenceType,Matches,{},\n' \
                        .format(','.join(sorted(sample[self.analysistype].allelenames)))
                    # Set the sequence counter to 0. This will be used when a sample has multiple best sequence types.
                    # The name of the sample will not be written on subsequent rows in order to make the report clearer
                    seqcount = 0
                    # Iterate through the best sequence types for the sample (only occurs if update profile is disabled)
                    for seqtype in self.resultprofile[sample.name]:
                        sample[self.analysistype].sequencetype = seqtype
                        # The number of matches to the profile
                        matches = self.resultprofile[sample.name][seqtype].keys()[0]
                        # If this is the first of one or more sequence types, include the sample name
                        if seqcount == 0:
                            try:
                                row += '{},{},{},{},'.format(sample.name, sample.mash.closestrefseqgenus, seqtype,
                                                             matches)
                            except KeyError:
                                row += '{},NA,{},{},'.format(sample.name, seqtype,
                                                             matches)
                        # Otherwise, skip the sample name
                        else:
                            row += ',,{},{},'.format(seqtype, matches)
                        # Iterate through all the genes present in the analyses for the sample
                        for gene in sorted(sample[self.analysistype].allelenames):
                            # refallele = self.profiledata[self.analysistype][seqtype][gene]
                            refallele = sample[self.analysistype].profiledata[seqtype][gene]
                            # Set the allele and percent id from the dictionary's keys and values, respectively
                            allele = self.resultprofile[sample.name][seqtype][matches][gene].keys()[0]
                            percentid = self.resultprofile[sample.name][seqtype][matches][gene].values()[0]
                            if refallele and refallele != allele:
                                if 0 < float(percentid) < 100:
                                    row += '{} ({}%) ({}),'.format(allele, percentid, refallele)
                                else:
                                    row += '{} ({}),'.format(allele, refallele)
                            else:
                                # Add the allele and % id to the row (only add the percent identity if it is not 100%)
                                if 0 < float(percentid) < 100:
                                    row += '{} ({}%),'.format(allele, percentid)
                                else:
                                    row += '{},'.format(allele)
                                    # self.referenceprofile[sample.name][gene] = allele
                        # Add a newline
                        row += '\n'
                        # Increment the number of sequence types observed for the sample
                        seqcount += 1
                    combinedrow += row
                    # If the length of the # of report directories is greater than 1 (script is being run as part of
                    # the assembly pipeline) make a report for each sample
                    with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                   self.analysistype), 'wb') as report:
                        # Write the row to the report
                        report.write(row)
                dotter()
            # Create the report folder
            make_path(self.reportpath)
            # Create the report containing all the data from all samples
            with open('{}/{}.csv'.format(self.reportpath, self.analysistype), 'wb') \
                    as combinedreport:
                # Write the results to this report
                combinedreport.write(combinedrow)
        # Clean up unnecessary metadata attributes
        self.metadatacleaner()

    def metadatacleaner(self):
        """Remove the attributes from the object; they take up too much room on the .json report"""
        for sample in self.runmetadata:
            try:
                delattr(sample[self.analysistype], "allelenames")
                delattr(sample[self.analysistype], "alleles")
                delattr(sample[self.analysistype], "profiledata")
                delattr(sample[self.analysistype], "outtables")
                delattr(sample[self.analysistype], "mashalleles")
                delattr(sample[self.analysistype], "faidict")
            except KeyError:
                pass
        printtime('{} analyses complete'.format(self.analysistype), self.start)

    def __init__(self, inputobject, analysistype):
        self.referencefilepath = inputobject.targetpath
        self.updatermlst = False
        self.alleles = ''
        self.profile = ''
        self.combinedalleles = ''
        self.supplementalprofile = ''
        self.profileset = set()
        self.cutoff = 0.8
        self.matchbonus = 2
        self.typequeue = Queue()
        self.allelefolders = set()
        self.bestmatch = defaultdict(int)
        self.bestdict = defaultdict(make_dict)
        self.mlstseqtype = defaultdict(make_dict)
        self.resultprofile = defaultdict(make_dict)
        self.targetpath = inputobject.targetpath
        self.customtargetpath = inputobject.customtargetpath
        self.reportpath = os.path.join(inputobject.path, 'reports')
        self.cpus = inputobject.cpus
        self.analysistype = analysistype
        self.sketchqueue = Queue(maxsize=self.cpus)
        self.mashqueue = Queue(maxsize=self.cpus)
        self.mashparsequeue = Queue(maxsize=self.cpus)
        self.databasesqueue = Queue(maxsize=self.cpus)
        self.baitqueue = Queue(maxsize=self.cpus)
        self.start = inputobject.starttime
        self.devnull = open(os.devnull, 'wb')
        # Custom.__init__(self, inputobject, analysistype, self.cutoff, self.matchbonus)
        self.runmetadata = inputobject.runmetadata
        # self.runmetadatacleaner()

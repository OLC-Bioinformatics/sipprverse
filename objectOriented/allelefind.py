#!/usr/bin/env python
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from customtargets import *
from glob import glob
__author__ = 'adamkoziol'


class FindAllele(Custom):

    def allelefinder(self):
        """
        Create an object for each strain to store metadata and get the name and the path of the reference allele
        files into the metadata object
        """
        printtime('Associating allele files with genes', self.start)
        for sample in self.metadata:
            # Create the desired attribute
            setattr(sample, self.analysistype, GenObject())
            # Initialise a dictionary to store name of the vtx gene and the path of the allele files
            sample[self.analysistype].targetfiles = dict()
            sample[self.analysistype].targets = dict()
            # For each vtx gene in the set, find the corresponding allele file
            for vtx in sample.general.vtxset:
                alleles = list()
                # Only the first four characters of the vtx gene are necessary when finding the corresponding allele
                # file e.g. stx1c -> stx1 will find the allele files for stx1A and stx1
                allelefiles = glob('{}/{}*.fa'.format(self.allelepath, vtx[:4]))
                for allelefile in allelefiles:
                    allele = os.path.basename(allelefile).split('.')[0]
                    # Populate the dictionaries
                    sample[self.analysistype].targetfiles[allele] = allelefile
                    alleles.append(allele)
                sample[self.analysistype].targets[vtx] = alleles
            sample[self.analysistype].targetpath = self.allelepath
        # Run the target baiting method
        self.targets()

    def targets(self):
        """Hash the allele files as required"""
        printtime('Performing analysis with {} targets folder'.format(self.analysistype), self.start)
        for sample in self.metadata:
            # Initialise dictionaries
            sample[self.analysistype].vtxhashes = dict()
            sample[self.analysistype].hashcalls = dict()
            try:
                # Iterate through all the vtx genes found for each strain
                for vtx, allelefile in sample[self.analysistype].targetfiles.items():
                    # Find the base name/path of the allele file
                    targetbase = allelefile.split('.')[0]
                    hashfile = '{}.mhs.gz'.format(targetbase)
                    # Define the hash call
                    hashcall = 'cd {} && mirabait -b {} -k 19 -K {}'.format(self.allelepath, allelefile, hashfile)
                    # Add the hash call to the dictionary
                    sample[self.analysistype].hashcalls[vtx] = hashcall
                    # Run the system call as required
                    if not os.path.isfile(hashfile):
                        call(hashcall, shell=True, stdout=self.devnull, stderr=self.devnull)
                    # Ensure that the hash file was successfully created
                    assert os.path.isfile(
                        hashfile), u'Hashfile could not be created for the combined target file {0!r:s}' \
                        .format(allelefile)
                    # Add the hash filename/path to the dictionary
                    sample[self.analysistype].vtxhashes[vtx] = hashfile
            except KeyError:
                pass
        # Bait the fastq files
        self.mirabaiting()

    def mirabaiting(self):
        """Perform baiting of fastq files with the appropriate hashed allele file"""
        printtime('Performing kmer baiting of fastq files with {} targets'.format(self.analysistype), self.start)
        # Create and start threads for each fasta file in the list
        for i in range(self.cpus):
            # Send the threads to the bait method
            threads = Thread(target=self.mirabait, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            try:
                # Create dictionaries to store data for each strain/vtx gene combination
                sample[self.analysistype].outputdir = dict()
                sample[self.analysistype].baitedfastq = dict()
                sample[self.analysistype].baitcall = dict()
                for vtx, vtxhash in sample[self.analysistype].vtxhashes.items():
                    # Set the output directory for each vtx gene
                    sample[self.analysistype].outputdir[vtx] = os.path.join(sample.general.outputdirectory,
                                                                            self.analysistype, vtx)
                    # Set attribute values
                    sample[self.analysistype].baitedfastq[vtx] = \
                        '{}/{}_targetMatches.fastq'.format(sample[self.analysistype].outputdir[vtx], vtx)
                    # Create the folder (if necessary)
                    make_path(sample[self.analysistype].outputdir[vtx])
                    # Create the mirabait call
                    if len(sample.general.fastqfiles) == 2:
                        syscall = 'mirabait -c -B {} -t 4 -o {} -p {} {}'\
                            .format(vtxhash, sample[self.analysistype].baitedfastq[vtx],
                                    sample.general.fastqfiles[0], sample.general.fastqfiles[1])
                    else:
                        syscall = 'mirabait -c -B {} -t 4 -o {} {}' \
                            .format(vtxhash, sample[self.analysistype].baitedfastq[vtx],
                                    sample.general.fastqfiles[0])
                    sample[self.analysistype].baitcall[vtx] = syscall
                    # Add the variables to the queue
                    self.baitqueue.put((sample, vtx))
            except KeyError:
                pass
        self.baitqueue.join()
        # Perform reference mapping
        self.mapping()

    def mirabait(self):
        while True:
            sample, vtx = self.baitqueue.get()
            # Run the system call (if necessary)
            if not os.path.isfile(sample[self.analysistype].baitedfastq[vtx]):
                call(sample[self.analysistype].baitcall[vtx], shell=True, stdout=self.devnull, stderr=self.devnull)
            self.baitqueue.task_done()

    def mapping(self):
        """Perform target indexing, reference mapping, SAM header editing, and BAM sorting using bowtie2 and samtools"""
        printtime('Performing reference mapping', self.start)
        for i in range(self.cpus):
            # Send the threads to
            threads = Thread(target=self.map, args=())
            # Set the daemon to True - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            try:
                # Create dictionaries to store data for each strain/vtx gene combination
                sample[self.analysistype].sortedbam = dict()
                sample[self.analysistype].baitfilenoext = dict()
                sample[self.analysistype].faifile = dict()
                sample[self.analysistype].bowtie2align = dict()
                sample[self.analysistype].bowtie2build = dict()
                sample[self.analysistype].samindex = dict()
                for vtx, vtxfastq in sample[self.analysistype].baitedfastq.items():
                    # Set the path/name for the sorted bam file to be created
                    sample[self.analysistype].sortedbam[vtx] = \
                        '{}/{}_sorted.bam'.format(sample[self.analysistype].outputdir[vtx], vtx)
                    # Remove the file extension of the bait file for use in the indexing command
                    sample[self.analysistype].baitfilenoext[vtx] = \
                        sample[self.analysistype].targetfiles[vtx].split('.')[0]
                    # Use bowtie2 wrapper to create index the target file
                    bowtie2build = Bowtie2BuildCommandLine(reference=sample[self.analysistype].targetfiles[vtx],
                                                           bt2=sample[self.analysistype].baitfilenoext[vtx],
                                                           **self.builddict)
                    # Use samtools wrapper to set up the bam sorting command
                    samsort = SamtoolsSortCommandline(input_bam=sample[self.analysistype].sortedbam[vtx],
                                                      o=True,
                                                      out_prefix="-")
                    # Create a list of programs to which data are piped as part of the reference mapping
                    samtools = [
                        # When bowtie2 maps reads to all possible locations rather than choosing a 'best' placement, the
                        # SAM header for that read is set to 'secondary alignment', or 256. Please see:
                        # http://davetang.org/muse/2014/03/06/understanding-bam-flags/ The script below reads stdin
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
                              '-U': sample[self.analysistype].baitedfastq[vtx],
                              '-a': True,
                              '--threads': self.cpus,
                              '--local': True}
                    # Create the bowtie2 reference mapping command
                    bowtie2align = Bowtie2CommandLine(bt2=sample[self.analysistype].baitfilenoext[vtx],
                                                      threads=self.cpus,
                                                      samtools=samtools,
                                                      **indict)
                    # Create the command to faidx index the bait file
                    sample[self.analysistype].faifile[vtx] = sample[self.analysistype].targetfiles[vtx] + '.fai'
                    samindex = SamtoolsFaidxCommandline(reference=sample[self.analysistype].targetfiles[vtx])
                    # Add the commands (as strings) to the metadata
                    sample[self.analysistype].bowtie2align[vtx] = str(bowtie2align)
                    sample[self.analysistype].bowtie2build[vtx] = str(bowtie2build)
                    sample[self.analysistype].samindex[vtx] = str(samindex)
                    # Index the allele files (if necessary)
                    if not os.path.isfile('{}.1{}'.format(sample[self.analysistype].baitfilenoext[vtx],
                                                          self.bowtiebuildextension)):
                        stdoutbowtieindex, stderrbowtieindex = \
                            map(StringIO, bowtie2build(cwd=sample[self.analysistype].targetpath))
                        # Write any error to a log file
                        if stderrbowtieindex:
                            # Write the standard error to log, bowtie2 puts alignment summary here
                            with open(os.path.join(sample[self.analysistype].targetpath,
                                                   '{}_bowtie_index.log'.format(vtx)), 'ab+') as log:
                                log.writelines(logstr(bowtie2build, stderrbowtieindex.getvalue(),
                                                      stdoutbowtieindex.getvalue()))
                        # Close the stdout and stderr streams
                        stdoutbowtieindex.close()
                        stderrbowtieindex.close()
                    # Add the commands to the queue. Note that the commands would usually be attributes of the sample
                    # but there was an issue with their serialization when printing out the metadata
                    self.mapqueue.put((sample, bowtie2build, bowtie2align, samindex, vtx))
            except KeyError:
                pass
        self.mapqueue.join()
        # Use samtools to index the sorted bam file
        self.indexing()

    def map(self):
        while True:
            sample, bowtie2build, bowtie2align, samindex, vtx = self.mapqueue.get()
            # Use samtools faidx to index the bait file - this will be used in the sample parsing
            if not os.path.isfile(sample[self.analysistype].faifile[vtx]):
                stdoutindex, stderrindex = map(StringIO, samindex(cwd=sample[self.analysistype].targetpath))
                # Write any error to a log file
                if stderrindex:
                    # Write the standard error to log, bowtie2 puts alignment summary here
                    with open(os.path.join(sample[self.analysistype].targetpath,
                                           '{}_samtools_index.log'.format(vtx)), 'ab+') as log:
                        log.writelines(logstr(samindex, stderrindex.getvalue(), stdoutindex.getvalue()))
                # Close the stdout and stderr streams
                stdoutindex.close()
                stderrindex.close()
            # Only run the functions if the sorted bam files and the indexed bait file do not exist
            if not os.path.isfile(sample[self.analysistype].sortedbam[vtx]):
                # Set stdout to a stringIO stream
                stdout, stderr = map(StringIO, bowtie2align(cwd=sample[self.analysistype].outputdir[vtx]))
                if stderr:
                    # Write the standard error to log, bowtie2 puts alignment summary here
                    with open(os.path.join(sample[self.analysistype].outputdir[vtx],
                                           '{}_bowtie_samtools.log'.format(vtx)), 'ab+') as log:
                        log.writelines(logstr(bowtie2align, stderr.getvalue(), stdout.getvalue()))
                stdout.close()
                stderr.close()
            self.mapqueue.task_done()

    def indexing(self):
        """Use samtools to index BAM files"""
        printtime('Indexing sorted bam files', self.start)
        for i in range(self.cpus):
            # Send the threads to
            threads = Thread(target=self.index, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            try:
                # Create dictionaries to store data for each strain/vtx gene combination
                sample[self.analysistype].bamindex = dict()
                sample[self.analysistype].sortedbai = dict()
                for vtx, sortedbam in sample[self.analysistype].sortedbam.items():
                    # Define the indexing call
                    bamindex = SamtoolsIndexCommandline(input=sortedbam)
                    # Update the metadata
                    sample[self.analysistype].sortedbai[vtx] = sortedbam + '.bai'
                    sample[self.analysistype].bamindex[vtx] = str(bamindex)
                    self.indexqueue.put((sample, bamindex, vtx))
            except KeyError:
                pass
        self.indexqueue.join()
        # Parse the results
        self.parsing()

    def index(self):
        while True:
            sample, bamindex, vtx = self.indexqueue.get()
            # Only make the call if the .bai file doesn't already exist
            if not os.path.isfile(sample[self.analysistype].sortedbai[vtx]):
                # Use cStringIO streams to handle bowtie output
                stdout, stderr = map(StringIO, bamindex(cwd=sample[self.analysistype].outputdir[vtx]))
                if stderr:
                    # Write the standard error to log
                    with open(os.path.join(sample[self.analysistype].outputdir[vtx],
                                           '{}_samtools_bam_index.log'.format(vtx)), 'ab+') as log:
                        log.writelines(logstr(bamindex, stderr.getvalue(), stdout.getvalue()))
                stderr.close()
            self.indexqueue.task_done()

    def parsing(self):
        printtime('Parsing sorted bam files', self.start)
        for i in range(self.cpus):
            # Send the threads to
            threads = Thread(target=self.parse, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Initialise dictionary-containing attributes that can be updated
            sample[self.analysistype].newsequences = dict()
            sample[self.analysistype].newseqclosestmatch = dict()
            sample[self.analysistype].newseqclosestseq = dict()
            sample[self.analysistype].allelematches = dict()
            try:
                for vtx, sortedbam in sample[self.analysistype].sortedbam.items():
                    self.parsequeue.put((sample, vtx))
            except KeyError:
                pass
        self.parsequeue.join()
        # Process the new alleles
        self.newalleles()

    def parse(self):
        import pysamstats
        import operator
        while True:
            sample, vtx = self.parsequeue.get()
            # Initialise dictionaries to store parsed data
            matchdict = dict()
            depthdict = dict()
            seqdict = dict()
            resultsdict = dict()
            snpdict = dict()
            gapdict = dict()
            faidict = dict()
            uniqueresults = dict()
            refdict = dict()
            # Variable to store the expected position in gene/allele
            pos = 0
            # Get the fai file into a dictionary to be used in parsing results
            with open(sample[self.analysistype].faifile[vtx], 'rb') as faifile:
                for line in faifile:
                    data = line.split('\t')
                    faidict[data[0]] = int(data[1])
            try:
                # Use the stat_variation function of pysam stats to return records parsed from sorted bam files
                # Values of interest can be retrieved using the appropriate keys
                correction = 0
                for rec in pysamstats.stat_variation(alignmentfile=sample[self.analysistype].sortedbam[vtx],
                                                     fafile=sample[self.analysistype].targetfiles[vtx],
                                                     max_depth=1000000):

                    # Add the reference sequence to the dictionary
                    if rec['chrom'] not in refdict:
                        refdict[rec['chrom']] = str()
                    refdict[rec['chrom']] += rec['ref']
                    # Initialise seqdict with the current gene/allele if necessary with an empty string
                    if rec['chrom'] not in seqdict:
                        seqdict[rec['chrom']] = str()
                        # Since this is the first position in a "new" gene/allele, reset the pos variable to 0
                        pos = 0
                        # There seems to be a bug in pysamstats with how gaps at the start of the sequence are treated.
                        # Although the position is correct, the whole reference sequence is still included, rather than
                        # starting at where the gap ends
                        if rec['pos'] > pos:
                            # If there is a gap of 173 bases at the beginning of the match, the reference sequence
                            # still should start at 0, but it starts at 173, therefore, the match actually starts at
                            # 2 * 173 = 346
                            correction = 2 * rec['pos']
                            # The number of gaps is equal to the starting position
                            gapdict[rec['chrom']] = rec['pos']
                            # The actual position will be rec['pos']
                            pos = rec['pos']
                    # Allow the position to reach the calculated correction factor
                    if rec['pos'] >= correction:
                        # Initialise gap dict with 0 gaps
                        if rec['chrom'] not in gapdict:
                            gapdict[rec['chrom']] = 0
                        # If there is a gap in the alignment, record the size of the gap in gapdict
                        if int(rec['pos']) > pos:
                            # Add the gap size to gap dict
                            gapdict[rec['chrom']] += rec['pos'] - pos
                            # Add dashes to the sequence to indicate the gap
                            seqdict[rec['chrom']] += 'N' * (int(rec['pos'] - pos))
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
                        # Track any deletions prior to the sequence
                        if rec['deletions'] > rec['matches']:
                            seqdict[rec['chrom']] += 'N'
                            # Increment the running total of the number of SNPs
                            snpdict[rec['chrom']] += 1
                        else:
                            if rec['matches'] > 0 or rec['mismatches'] > 0:
                                # If the most prevalent base (calculated with max() and operator.itemgetter())
                                # doesn't match the reference base, add this prevalent base to seqdict
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
                # Calculate the average depth by dividing the total number of reads observed by the length of the gene
                averagedepth = float(depthdict[allele]) / float(matchdict[allele])
                percentidentity = float(matchdict[allele]) / float(faidict[allele]) * 100
                if percentidentity == self.cutoff * 100:
                    # Only report a positive result if this average depth is greater than 4X
                    if averagedepth > 4:
                        # Populate resultsdict with the gene/allele name, the percent identity, and the average depth
                        resultsdict.update({allele: {'{:.2f}'.format(percentidentity): '{:.2f}'.format(averagedepth)}})
            # Add the results to the object
            sample[self.analysistype].allelematches[vtx] = resultsdict
            # Determine if there are alleles without a 100% match
            if not resultsdict:
                for allele in sorted(matchdict):
                    # Filter the alleles to only include the vtx subunit
                    if vtx in allele:
                        percentidentity = float(matchdict[allele]) / float(faidict[allele]) * 100
                        # Use a more relaxed cutoff to find the closest alleles
                        if percentidentity >= self.cutoff * 50:
                            uniqueresults.update({allele: percentidentity})
                try:
                    # Find the best match (highest percent identity)
                    closestallele = max(uniqueresults.iteritems(), key=operator.itemgetter(1))[0]
                    percentidentity = max(uniqueresults.iteritems(), key=operator.itemgetter(1))[1]
                    averagedepth = float(depthdict[closestallele]) / float(matchdict[closestallele])
                    # Populate the metadata with the results
                    sample[self.analysistype].newsequences[vtx] = seqdict[closestallele]
                    sample[self.analysistype].newseqclosestmatch[vtx] = \
                        {closestallele: {'{:.2f}'.format(percentidentity): '{:.2f}'.format(averagedepth)}}
                    sample[self.analysistype].newseqclosestseq[vtx] = {closestallele: refdict[closestallele]}
                except ValueError:
                    pass
            self.parsequeue.task_done()

    def newalleles(self):
        """Determine whether 'new' alleles have been previously stored in the custom database"""
        from Bio import SeqIO
        from Bio.Alphabet import generic_dna
        printtime('Updating metadata with new allele information', self.start)
        # Iterate through all the samples
        for sample in self.metadata:
            # Create a dictionary to store brand new alleles
            sample[self.analysistype].allelesnew = dict()
            sample[self.analysistype].allelesolc = dict()
            # Find the vtx gene and sequence for new alleles
            for vtx, sequence in sample[self.analysistype].newsequences.items():
                # Extract the subtype for the vtx gene from the metadata
                subtype = [stx[-1] for stx in sample.general.vtxset if vtx[:4] == stx[:4]][0]
                # Populate a dictionary with the sequence as the key - duplicates will automatically be discarded
                self.alleledict[sequence] = {vtx: subtype}
                # Set the file to in which the new allele is to be stored
                allelefile = '{}/OLC{}.tfa'.format(sample[self.analysistype].targetpath, vtx)
                # Create the file if it doesn't exist
                open(allelefile, 'ab').close()
                self.allelefiledict[vtx] = allelefile
        # Iterate through all the unique new alleles to find alleles not already in the custom database
        for sequence, vtxtype in self.alleledict.items():
            vtx, subtype = vtxtype.items()[0]
            # Load the allele file into a list
            allelelist = list(SeqIO.parse(self.allelefiledict[vtx], 'fasta'))
            # Put the most recent allele number (+1) for each gene into a dictionary
            try:
                self.allelelist.update({vtx: int(allelelist[-1].id.split(':')[0]) + 1})
            # If there are no previous entries, start at 1000000
            except IndexError:
                self.allelelist.update({vtx: 1000000})
            # List comprehension to store all the sequence of all the alleles
            alleles = [allele.seq for allele in allelelist]
            # If this sequence has not previously been found add it to the dictionary
            if sequence not in alleles:
                self.newallelelist[sequence] = {vtx: subtype}
            # If it has been found before, update the metadata to reflect this
            else:
                for sample in self.metadata:
                    currentid = allelelist[-1].id
                    # The .newsequences attribute is to be updated
                    for stx, seq in sample[self.analysistype].newsequences.items():
                        # If the current sequence is the same as the sequence store in the metadata object
                        if sequence == seq:
                            # Pull the depth of coverage value from the .newseqclosestmatch attribute
                            depth = sample[self.analysistype].newseqclosestmatch[vtx].items()[0][1].items()[0][1]
                            # Update the metadata
                            sample[self.analysistype].allelematches[stx] = {currentid: {100.00: depth}}
                            # Add metadata with any custom alleles
                            sample[self.analysistype].allelesolc[vtx] = {'nt': {currentid: sequence}}
                            # Perform translations to quickly find any alleles with mis/nonsense mutations
                            # Add an appropriate number of N's to pad out any partial codons at the end of the sequence
                            remainder = 3 - len(sequence) % 3
                            sequence += ('N' * remainder)
                            protein = Seq.translate(Seq(sequence))
                            sample[self.analysistype].allelesolc[vtx].update({'aa': {currentid: str(protein)}})
        # Add the alleles to the appropriate file
        for sequence, vtxtype in self.newallelelist.items():
            # Split the tuple
            vtx, subtype = vtxtype.items()[0]
            try:
                # The header will be the allele number plus the predicted subtype
                currentid = '{}:{}'.format(str(self.allelelist[vtx]), subtype)
                allelesequence = Seq(sequence, generic_dna)
                # Create a sequence record using BioPython
                fasta = SeqRecord(allelesequence,
                                  # Without this, the header will be improperly formatted
                                  description='',
                                  # Use >:currentid as the header
                                  id=currentid)
                # Open the allele file to append
                with open(self.allelefiledict[vtx], 'ab+') as supplemental:
                    # Use the SeqIO module to properly format the new sequence record
                    SeqIO.write(fasta, supplemental, "fasta")
                # Populate the metadata with the new allele information
                for sample in self.metadata:
                    for stx, seq in sample[self.analysistype].newsequences.items():
                        if sequence == seq:
                            # Pull the depth of coverage value from the .newseqclosestmatch attribute
                            depth = sample[self.analysistype].newseqclosestmatch[vtx].items()[0][1].items()[0][1]
                            # Update the metadata
                            sample[self.analysistype].allelematches[stx] = {currentid: {100.00: depth}}
                            sample[self.analysistype].allelesnew[vtx] = currentid
                            # Add metadata with any custom alleles
                            sample[self.analysistype].allelesolc[vtx] = {'nt': {currentid: sequence}}
                            # Perform translations to quickly find any alleles with mis/nonsense mutations
                            # Add an appropriate number of N's to pad out any partial codons at the end of the sequence
                            remainder = 3 - len(sequence) % 3
                            sequence += ('N' * remainder)
                            protein = Seq.translate(Seq(sequence))
                            sample[self.analysistype].allelesolc[vtx].update({'aa': {currentid: str(protein)}})
                # Increment the current id
                self.allelelist[vtx] += 1
            except KeyError:
                pass
        # Create a report of all the new alleles
        self.reports()

    def reports(self):
        """Create a report with all the new alleles"""
        import xlsxwriter
        printtime('Creating report', self.start)
        # If, for some reason, analyses are performed more than once on this dataset, alleles will no longer be
        # considered 'new' on subsequent analyses. Don't create/overwrite new allele reports if there are no new alleles
        if self.newallelelist:
            # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
            # able to have appropriately sized, multi-line cells
            workbook = xlsxwriter.Workbook('{}/newalleles.xlsx'.format(self.reportpath))
            # New worksheet to store the data
            worksheet = workbook.add_worksheet()
            # Add a bold format for header cells. Using a monotype font size 8
            bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 8})
            # Format for data cells. Monotype, size 8, top vertically justified
            courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 8})
            courier.set_align('top')
            # Set the custom width for columns 3 and 4 to be 50, and column 5 to be 60 characters of the default font
            worksheet.set_column(3, 4, 50)
            worksheet.set_column(5, 5, 60)
            # Initialise the position within the worksheet to be (0,0)
            row = 0
            col = 0
            # List of the headers to use
            headers = ['Strain', 'Subunit', 'Closest', 'NucleicAcidSequence', 'ProteinSequence', 'BLASTalignment']
            # Populate the headers
            for category in headers:
                # Write the data in the specified cell (row, col) using the bold format
                worksheet.write(row, col, category, bold)
                # Move to the next column to write the next category
                col += 1
            # Data starts in row 1
            row = 1
            # Initialise variables to hold the longest names; used in setting the column width
            longestname = 0
            longestrefname = 0
            for sample in self.metadata:
                if sample[self.analysistype].allelesnew.items():
                    # Every record starts at column 0
                    col = 0
                    # List to store data; will be used for populating the spreadsheet
                    results = list()
                    for vtx, allele in sample[self.analysistype].allelesnew.items():
                        # Extract the protein sequence of the gene
                        protein = sample[self.analysistype].allelesolc[vtx]['aa'][allele]
                        # Format the protein as fasta
                        aastring = self.fastacarriagereturn(allele, protein, '\n')
                        # Extract the nucleotide sequence, and format it to fasta
                        ntsequence = sample[self.analysistype].allelesolc[vtx]['nt'][allele]
                        ntstring = self.fastacarriagereturn(allele, ntsequence, '\n')
                        # Write the new alleles to a fasta file
                        fastafile = '{}/newalleles.fa'.format(self.reportpath)
                        with open(fastafile, 'wb') as newalleles:
                            newalleles.write(ntstring)
                        # Extract the name of the closest reference allele, and its sequence
                        refallele, refsequence = sample[self.analysistype].newseqclosestseq[vtx].items()[0]
                        # Get the percent identity of this closest reference allele
                        percentidentity = sample[self.analysistype].newseqclosestmatch[vtx].items()[0][1].items()[0][0]
                        # Create a pseudo-BLAST alignment of the query and reference sequences
                        formattedblast = self.interleaveblastresults(ntsequence, refsequence)
                        # Determine the longest name of all the strains, and use it to set the width of column 0
                        if len(sample.name) > longestname:
                            longestname = len(sample.name)
                            worksheet.set_column(0, 0, len(sample.name))
                        # Do the same for the reference allele names
                        if len(refallele) > longestrefname:
                            longestrefname = len(refallele)
                            worksheet.set_column(2, 2, longestrefname)
                        # Set the width of the row to be the number of lines (number of newline characters) * 11
                        worksheet.set_row(row, formattedblast.count('\n') * 11)
                        # Store the variables in a list
                        results = [sample.name, vtx, '{}\n\n{}%'.format(refallele, percentidentity), ntstring, aastring,
                                   formattedblast]
                    # Write out the data to the spreadsheet
                    for data in results:
                        worksheet.write(row, col, data, courier)
                        col += 1
                    # Increase the row counter for the next strain's data
                    row += 1
            # Close the workbook
            workbook.close()

    @staticmethod
    def interleaveblastresults(query, subject):
        """
        Creates an interleaved string that resembles BLAST sequence comparisons
        :param query: Query sequence
        :param subject: Subject sequence
        :return: Properly formatted BLAST-like sequence comparison
        """
        # Initialise strings to hold the matches, and the final BLAST-formatted string
        matchstring = ''
        blaststring = ''
        # Iterate through the query
        for i, bp in enumerate(query):
            # If the current base in the query is identical to the corresponding base in the reference, append a '|'
            # to the match string, otherwise, append a ' '
            if bp == subject[i]:
                matchstring += '|'
            else:
                matchstring += ' '
        # Set a variable to store the progress through the sequence
        prev = 0
        # Iterate through the query, from start to finish in steps of 60 bp
        for j in range(0, len(query), 60):
            # BLAST results string. The components are: current position (padded to four characters), 'OLC', query
            # sequence, \n, matches, \n, 'ref', subject sequence. Repeated until all the sequence data are present.
            """
            0000 OLC ATGAAGAAGATATTTGTAGCGGCTTTATTTGCTTTTGTTTCTGTTAATGCAATGGCAGCT
                     ||||||||||| ||| | |||| ||||||||| || ||||||||||||||||||||||||
                 ref ATGAAGAAGATGTTTATGGCGGTTTTATTTGCATTAGTTTCTGTTAATGCAATGGCAGCT
            0060 OLC GATTGTGCAAAAGGTAAAATTGAGTTCTCTAAGTATAATGAGAATGATACATTCACAGTA
                     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                 ref GATTGTGCAAAAGGTAAAATTGAGTTCTCTAAGTATAATGAGAATGATACATTCACAGTA
            """
            blaststring += '{} OLC {}\n         {}\n     ref {}\n'\
                .format('{:04d}'.format(j), query[prev:j + 60], matchstring[prev:j + 60], subject[prev:j + 60])
            # Update the progress variable
            prev = j + 60
        # Return the properly formatted string
        return blaststring

    @staticmethod
    def fastacarriagereturn(seqname, seq, delim='\015'):
        """
        Format sequences to be fasta with carriage returns ('\015') instead of newlines
        :param seqname: Name of sequence
        :param seq: Sequence string
        :param delim: Delimiter used to split fasta sequence. Defaults to carriage return
        :return: string of properly formatted sequence
        """
        count = 0
        # Standard fasta header
        seqstring = '>{}{}'.format(seqname, delim)
        # Iterate through the sequence, adding a delimiter every 60 characters
        for char in seq:
            if count < 60:
                seqstring += char
                count += 1
            else:
                seqstring += char + delim
                count = 0
        # Add a final delimiter at the end of the string
        seqstring += delim
        # Return the properly formatted string
        return seqstring

    def __init__(self, inputobject, analysistype):
        # Create a custom object using a cutoff of 100%
        Custom.__init__(self, inputobject, analysistype, 1.00)
        # Initialise variables
        self.allelepath = os.path.join(self.targetpath, 'vtxalleles')
        self.allelefiles = glob('{}/*.fa'.format(self.allelepath))
        self.alleledict = dict()
        self.allelefiledict = dict()
        self.allelelist = dict()
        self.newallelelist = dict()
        self.reportpath = os.path.join(self.path, 'reports')

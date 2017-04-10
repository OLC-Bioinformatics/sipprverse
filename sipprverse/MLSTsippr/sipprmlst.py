#!/usr/bin/env python
from sipprcommon.sippingmethods import *
__author__ = 'adamkoziol'


class MLSTmap(Sippr):

    def targets(self):
        from Bio import SeqIO
        printtime('Finding {} target files'.format(self.analysistype), self.start)
        for sample in self.runmetadata:
            setattr(sample, self.analysistype, GenObject())
            if not self.pipeline:
                self.profile = glob('{}*.txt'.format(self.targetpath))
                self.combinedalleles = glob('{}/*.fasta'.format(self.targetpath))[0]
                geneset = set()
                if self.analysistype.lower() == 'rmlst':
                    # Find all the gene names from the combined alleles files
                    for record in SeqIO.parse(open(self.combinedalleles, "rU"), "fasta"):
                        # Add them to the set of names
                        geneset.add(record.id.split('_')[0])
                else:
                    # Find all the gene names from the combined alleles files
                    for record in SeqIO.parse(open(self.combinedalleles, "rU"), "fasta"):
                        # Add them to the set of names
                        geneset.add(record.id.split('-')[0])
                # Add the combined alleles to the profile set
                self.profileset.add(self.combinedalleles)
                sample[self.analysistype].targetpath = self.targetpath
                sample[self.analysistype].alleledir = self.targetpath
                sample[self.analysistype].alleles = sorted(list(geneset))
                sample[self.analysistype].allelenames = sorted([os.path.split(x)[1].split('.')[0] for x in
                                                                sample[self.analysistype].alleles])
                sample[self.analysistype].profile = self.profile
                sample[self.analysistype].analysistype = self.analysistype
                sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory, self.analysistype)
                sample[self.analysistype].combinedalleles = self.combinedalleles
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
            # Populate the appropriate attributes
            for sample in self.runmetadata:
                if sample.general.bestassemblyfile != 'NA':
                    if sample[self.analysistype].combinedalleles == target:
                        sample[self.analysistype].hashcall = hashcall
                        sample[self.analysistype].hashfile = hashfile
                        sample[self.analysistype].targetpath = self.targetpath
                        sample[self.analysistype].outputdir = os.path.join(sample.run.outputdirectory, self.analysistype)
                        sample[self.analysistype].baitedfastq = '{}/{}_targetMatches.fastq'\
                            .format(sample[self.analysistype].outputdir, self.analysistype)
        # Run the baiting method in the Sippr class
        self.baiting()

    def __init__(self, inputobject, analysistype):
        self.analysistype = analysistype
        self.targetpath = inputobject.targetpath
        self.profile = list()
        self.combinedalleles = list()
        self.profileset = set()
        self.runmetadata = inputobject.runmetadata.samples
        self.pipeline = inputobject.pipeline
        self.copy = inputobject.copy
        Sippr.__init__(self, inputobject)

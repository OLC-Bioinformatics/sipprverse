#!/usr/bin/env python
from sipprCommon.sippingmethods import Sippr
from accessoryFunctions.accessoryFunctions import printtime, GenObject
import os
from glob import glob
__author__ = 'adamkoziol'


class MLSTmap(Sippr):

    def targets(self):
        from Bio import SeqIO
        printtime('Finding {} target files'.format(self.analysistype), self.start)
        for sample in self.runmetadata:
            setattr(sample, self.analysistype, GenObject())
            sample[self.analysistype].runanalysis = True
            if not self.pipeline:
                # self.profile = glob('{}*.txt'.format(self.targetpath))
                # self.combinedalleles = glob('{}/*.fasta'.format(self.targetpath))[0]
                try:
                    sample[self.analysistype].profile = glob(os.path.join(self.targetpath, '*.txt'))[0]
                    sample[self.analysistype].combinedalleles = glob(os.path.join(self.targetpath, '*.fasta'))[0]
                    sample[self.analysistype].targetpath = self.targetpath
                    sample[self.analysistype].alleledir = self.targetpath
                except IndexError:
                    print('Cannot find the profile and/or the combined allele file in the designated target path ({})'
                          'please ensure that those files exist'.format(self.targetpath))
                    quit()
            else:

                try:
                    if self.analysistype.lower() == 'rmlst':
                        targetdir = self.targetpath
                    else:
                        targetdir = os.path.join(self.targetpath, sample.general.closestrefseqgenus)
                    sample[self.analysistype].profile = glob(os.path.join(targetdir, '*.txt'))[0]
                    sample[self.analysistype].combinedalleles = glob(os.path.join(targetdir, '*.fasta'))[0]
                    sample[self.analysistype].targetpath = targetdir
                    sample[self.analysistype].alleledir = targetdir
                except IndexError:
                    sample[self.analysistype].profile = 'NA'
                    sample[self.analysistype].combinedalleles = 'NA'
                    sample[self.analysistype].runanalysis = False
            geneset = set()
            try:
                # Find all the gene names from the combined alleles files
                for record in SeqIO.parse(open(sample[self.analysistype].combinedalleles, "rU"), "fasta"):
                    # Determine whether an underscore, or a hyphen is being used to separate the gene name and allele
                    # number, split on the delimiter, and add the gene name to the set
                    if '_' in record.id:
                        geneset.add(record.id.split('_')[0])
                    elif '-' in record.id:
                        geneset.add(record.id.split('-')[0])
            except FileNotFoundError:
                sample[self.analysistype].runanalysis = False

            # Add the combined alleles to the profile set
            self.profileset.add(sample[self.analysistype].combinedalleles)
            sample[self.analysistype].alleles = sorted(list(geneset)) if geneset else 'NA'
            sample[self.analysistype].allelenames = \
                sorted([os.path.splitext(os.path.split(x)[1])[0] for x in sample[self.analysistype].alleles]) if \
                    geneset else 'NA'
            #
            sample[self.analysistype].analysistype = self.analysistype
            sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory, self.analysistype)
            sample[self.analysistype].baitfile = sample[self.analysistype].combinedalleles

        printtime('Indexing {} target file'.format(self.analysistype), self.start)
        # Ensure that the hash file was successfully created
        # Populate the appropriate attributes
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                # Set the necessary attributes
                sample[self.analysistype].outputdir = os.path.join(sample.general.outputdirectory,
                                                                   self.analysistype)
                sample[self.analysistype].logout = os.path.join(sample[self.analysistype].outputdir, 'logout.txt')
                sample[self.analysistype].logerr = os.path.join(sample[self.analysistype].outputdir, 'logerr.txt')
                sample[self.analysistype].baitedfastq = '{}/{}_targetMatches.fastq.gz'\
                    .format(sample[self.analysistype].outputdir, self.analysistype)
        # Run the baiting method in the Sippr class
        self.bait()

    def __init__(self, inputobject, analysistype, cutoff):
        self.analysistype = analysistype
        self.targetpath = inputobject.targetpath
        self.profileset = set()
        self.runmetadata = inputobject.runmetadata.samples
        self.pipeline = inputobject.pipeline
        self.copy = inputobject.copy
        self.logfile = inputobject.logfile
        Sippr.__init__(self, inputobject, cutoff)

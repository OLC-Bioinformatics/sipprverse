#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import printtime, make_path, GenObject
from threading import Thread
from subprocess import call
from queue import Queue
import os
import re
__author__ = 'adamkoziol'


class Mash(object):
    def sketching(self):
        printtime('Indexing files for {} analysis'.format(self.analysistype), self.starttime)
        # Create the threads for the analysis
        for i in range(self.cpus):
            threads = Thread(target=self.sketch, args=())
            threads.setDaemon(True)
            threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            # Create the analysis type-specific GenObject
            setattr(sample, self.analysistype, GenObject())
            # Set attributes
            sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory, self.analysistype)
            make_path(sample[self.analysistype].reportdir)
            sample[self.analysistype].targetpath = self.referencefilepath if not self.pipeline else os.path.join(
                self.referencefilepath, self.analysistype)
            sample[self.analysistype].refseqsketch = os.path.join(sample[self.analysistype].targetpath,
                                                                  'RefSeqSketchesDefaults.msh')
            sample[self.analysistype].sketchfilenoext = os.path.join(sample[self.analysistype].reportdir, sample.name)
            sample[self.analysistype].sketchfile = sample[self.analysistype].sketchfilenoext + '.msh'
            # Make the mash output directory if necessary
            make_path(sample[self.analysistype].reportdir)
            # Create a file containing the path/name of the filtered, corrected fastq files
            sample[self.analysistype].filelist = os.path.join(sample[self.analysistype].reportdir,
                                                              '{}_fastqfiles.txt'.format(sample.name))
            with open(sample[self.analysistype].filelist, 'w') as filelist:
                filelist.write('\n'.join(sample.general.trimmedcorrectedfastqfiles))

            # Create the system call
            sample.commands.sketch = 'mash sketch -m 2 -p {} -l {} -o {}' \
                .format(self.cpus, sample[self.analysistype].filelist, sample[self.analysistype].sketchfilenoext)
            # Add each sample to the threads
            try:
                self.sketchqueue.put(sample)
            except (KeyboardInterrupt, SystemExit):
                printtime('Received keyboard interrupt, quitting threads', self.starttime)
                quit()
        # Join the threads
        self.sketchqueue.join()
        self.mashing()

    def sketch(self):
        while True:
            sample = self.sketchqueue.get()
            if not os.path.isfile(sample[self.analysistype].sketchfile):
                call(sample.commands.sketch, shell=True, stdout=self.fnull, stderr=self.fnull)
            self.sketchqueue.task_done()

    def mashing(self):
        printtime('Performing {} analyses'.format(self.analysistype), self.starttime)
        # Create the threads for the analysis
        for i in range(self.cpus):
                threads = Thread(target=self.mash, args=())
                threads.setDaemon(True)
                threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            sample[self.analysistype].mashresults = os.path.join(sample[self.analysistype].reportdir, '{}.tab'.format(
                sample.name))

            sample.commands.mash = \
                'mash dist -p {} {} {} | sort -gk3 > {}'.format(self.threads,
                                                                sample[self.analysistype].refseqsketch,
                                                                sample[self.analysistype].sketchfile,
                                                                sample[self.analysistype].mashresults)
            try:
                self.mashqueue.put(sample)
            except (KeyboardInterrupt, SystemExit):
                printtime('Received keyboard interrupt, quitting threads', self.starttime)
                quit()
        # Join the threads
        self.mashqueue.join()
        self.parse()

    def mash(self):
        while True:
            sample = self.mashqueue.get()
            #
            if not os.path.isfile(sample[self.analysistype].mashresults):
                call(sample.commands.mash, shell=True, stdout=self.fnull, stderr=self.fnull)
            self.mashqueue.task_done()

    def parse(self):
        printtime('Determining closest refseq genome', self.starttime)
        # Create a dictionary to store the accession: taxonomy id of refseq genomes
        refdict = dict()
        # Set the name of the file storing the assembly summaries
        referencefile = os.path.join(self.referencefilepath, self.analysistype, 'assembly_summary_refseq.txt')
        # Extract the accession: genus species key: value pairs from the refseq summary file
        with open(referencefile) as reffile:
            for line in reffile:
                # Ignore the first couple of lines
                if line.startswith('# assembly_accession'):
                    # Iterate through all the lines with data
                    for accessionline in reffile:
                        # Replace commas with semicolons
                        accessionline = accessionline.replace(',', ';')
                        # Split the lines on tabs
                        data = accessionline.split('\t')
                        # Populate the dictionary with the accession: tax id e.g. GCF_001298055: Helicobacter pullorum
                        refdict[data[0].split('.')[0]] = data[7]
        for sample in self.metadata:
            # Initialise a list to store all the MASH results
            mashdata = list()
            # Open the results and extract the data
            with open(sample[self.analysistype].mashresults, 'r') as results:
                for line in results:
                    mashdata.append(line.rstrip())
            # Ensure that there is at least a single result
            if mashdata:
                # Iterate through the data
                for linedata in mashdata:
                    try:
                        # Split on tabs
                        data = linedata.split('\t')
                        # Extract the components of the line
                        referenceid, queryid, sample[self.analysistype].mashdistance, sample[self.analysistype]. \
                            pvalue, sample[self.analysistype].nummatches = data
                        # Extract the name of the refseq assembly from the mash outputs, and split as necessary e.g.
                        # GCF_000008865.1_ASM886v1_genomic.fna.gz becomes GCF_000008865
                        refid = referenceid.split('.')[0]
                        # Find the genus and species of the sample using the dictionary of refseq summaries
                        sample[self.analysistype].closestrefseq = refdict[refid]
                        sample[self.analysistype].closestrefseqgenus = \
                            sample[self.analysistype].closestrefseq.split()[0]
                        sample[self.analysistype].closestrefseqspecies = \
                            sample[self.analysistype].closestrefseq.split()[1]
                        # Set the closest refseq genus - will be used for all typing that requires the genus to be known
                        sample.general.referencegenus = sample[self.analysistype].closestrefseqgenus
                        break
                    except ValueError:
                        sample[self.analysistype].closestrefseq = 'NA'
                        sample[self.analysistype].closestrefseqgenus = 'NA'
                        sample[self.analysistype].closestrefseqspecies = 'NA'
                        sample[self.analysistype].mashdistance = 'NA'
                        sample[self.analysistype].pvalue = 'NA'
                        sample[self.analysistype].nummatches = 'NA'
                        sample.general.referencegenus = 'NA'
                        break
                    # I have encountered a strain that has a best match with the RefSeq database provided by the
                    # developers of MASH that doesn't have a corresponding entry in the assembly_summary_refseq.txt
                    # file downloaded from NCBI. Simply pass on this and look for the next best hit
                    except KeyError:
                        pass
            else:
                sample[self.analysistype].closestrefseq = 'NA'
                sample[self.analysistype].closestrefseqgenus = 'NA'
                sample[self.analysistype].closestrefseqspecies = 'NA'
                sample[self.analysistype].mashdistance = 'NA'
                sample[self.analysistype].pvalue = 'NA'
                sample[self.analysistype].nummatches = 'NA'
                sample.general.referencegenus = 'NA'
        self.reporter()

    def reporter(self):
        """
        Create the MASH report
        """
        printtime('Creating {} report'.format(self.analysistype), self.starttime)
        make_path(self.reportpath)
        header = 'Strain,ReferenceGenus,ReferenceFile,ReferenceGenomeMashDistance,Pvalue,NumMatchingHashes\n'
        data = ''
        for sample in self.metadata:
            try:
                data += '{},{},{},{},{},{}\n'.format(sample.name,
                                                     sample[self.analysistype].closestrefseqgenus,
                                                     sample[self.analysistype].closestrefseq,
                                                     sample[self.analysistype].mashdistance,
                                                     sample[self.analysistype].pvalue,
                                                     sample[self.analysistype].nummatches)
            except AttributeError:
                data += '{}\n'.format(sample.name)
        # Create the report file
        reportfile = os.path.join(self.reportpath, 'mash.csv')
        with open(reportfile, 'w') as report:
            report.write(header)
            report.write(data)

    def __init__(self, inputobject, analysistype):
        self.metadata = inputobject.runmetadata.samples
        self.referencefilepath = inputobject.reffilepath
        self.starttime = inputobject.starttime
        self.reportpath = inputobject.reportpath
        self.cpus = inputobject.cpus
        self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(self.metadata) > 1 else 1
        self.sketchqueue = Queue(maxsize=self.cpus)
        self.mashqueue = Queue(maxsize=4)
        self.analysistype = analysistype
        self.pipeline = inputobject.pipeline
        self.fnull = open(os.devnull, 'w')  # define /dev/null
        self.sketching()

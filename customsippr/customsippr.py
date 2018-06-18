#!/usr/bin/python3
from accessoryFunctions.accessoryFunctions import combinetargets, printtime
from sipprCommon.sippingmethods import Sippr
from Bio import SeqIO
import os

__author__ = 'adamkoziol'


class CustomGenes(object):

    def main(self):
        """
        Run the necessary methods in the correct order
        """
        self.target_validate()
        self.gene_names()
        Sippr(self)
        self.report()

    def target_validate(self):
        """
        Validate the user-supplied targets by running the (multi-)FASTA file through the method that combines the
        targets. Will also be useful for the downstream analyses
        """
        printtime('Validating user-supplied targets', self.starttime)
        combinetargets([self.targets], self.targetpath)

    def gene_names(self):
        """
        Extract the names of the user-supplied targets
        """
        # Iterate through all the target names in the formatted targets file
        for record in SeqIO.parse(self.targets, 'fasta'):
            # Append all the gene names to the list of names
            self.genes.append(record.id)

    def report(self):
        """
        Create the report for the user-supplied targets
        """
        # Add all the genes to the header
        header = 'Sample,{genes}\n'.format(genes=','.join(sorted(self.genes)))
        data = str()
        with open(os.path.join(self.reportpath, '{}.csv'.format(self.analysistype)), 'w') as report:
            for sample in self.runmetadata:
                data += sample.name + ','
                # Iterate through all the user-supplied target names
                for target in sorted(self.genes):
                    gene_results = str()
                    for name, identity in sample[self.analysistype].results.items():
                        # If the current target matches the target in the header, add the data to the string
                        if name == target:
                            gene_results += '{percent_id}% ({avgdepth} +/ {stddev}),'\
                                .format(percent_id=identity,
                                        avgdepth=sample[self.analysistype].avgdepth[name],
                                        stddev=sample[self.analysistype].standarddev[name])
                    # Populate the data string appropriately
                    if gene_results:
                        data += gene_results
                    else:
                        data += ','
                data += '\n'
            # Write the strings to the report
            report.write(header)
            report.write(data)

    def __init__(self, args):
        self.starttime = args.starttime
        self.targets = args.user_genes
        self.targetpath = os.path.split(self.targets)[0]
        self.path = args.path
        self.reportpath = args.reportpath
        self.runmetadata = args.runmetadata.samples
        self.start = args.starttime
        self.analysistype = 'custom'
        self.cpus = args.cpus
        self.threads = args.threads
        self.homepath = args.homepath
        self.sequencepath = args.sequencepath
        self.pipeline = False
        self.taxonomy = False
        self.logfile = args.logfile
        self.genes = list()

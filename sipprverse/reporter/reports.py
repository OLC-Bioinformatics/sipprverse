#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import make_path, MetadataObject, printtime
from glob import glob
import numpy
import os
__author__ = 'adamkoziol'


class Reports(object):

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

    def __init__(self, inputobject):
        self.starttime = inputobject.starttime
        self.samples = inputobject.samples
        self.completemetadata = inputobject.completemetadata
        self.path = inputobject.path
        self.analysescomplete = inputobject.analysescomplete
        self.reportpath = inputobject.reportpath
        self.runmetadata = inputobject.runmetadata
        self.analysistype = inputobject.analysistype

#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import make_path, MetadataObject, printtime
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from glob import glob
import operator
import numpy
import os
__author__ = 'adamkoziol'


class Reports(object):

    def reporter(self, analysistype='genesippr'):
        """
        Creates a report of the genesippr results
        :param analysistype: The variable to use when accessing attributes in the metadata object
        """
        printtime('Creating {} report'.format(analysistype), self.starttime, output=self.portallog)
        # Create a dictionary to link all the genera with their genes
        genusgenes = dict()
        # A list to store all the unique gene names
        geneset = list()
        # The organism-specific targets are in .tfa files in the target path
        targetpath = str()
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                targetpath = sample[analysistype].targetpath
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
            sample[analysistype].targetgenera = list()
            if sample.general.bestassemblyfile != 'NA':
                for organism in genusgenes:
                    # Iterate through all the genesippr hits and attribute each gene to the appropriate genus
                    for gene in sample[analysistype].results:
                        # If the gene name is in the genes from that organism, add the genus name to the list of
                        # genera found in the sample
                        if gene in genusgenes[organism]:
                            if organism not in sample[analysistype].targetgenera:
                                sample[analysistype].targetgenera.append(organism)
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # The report will have every gene for all genera in the header
        header = 'Strain,Genus,{},\n'.format(','.join(geneset))
        data = str()
        with open(os.path.join(self.reportpath, analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                if sample.general.bestassemblyfile != 'NA':
                    # Add the genus/genera found in the sample
                    data += '{},{},'.format(sample.name, ';'.join(sample[analysistype].targetgenera))
                    if sample[analysistype].results:
                        for gene in geneset:
                            # If the gene was not found in the sample, print an empty cell in the report
                            if gene not in [target[0] for target in sample[analysistype].results.items()]:
                                data += ','
                            # Print the required information for the gene
                            for name, identity in sample[analysistype].results.items():
                                if name == gene:
                                    data += '{}% ({} +/- {}),'.format(identity,
                                                                      sample[analysistype].avgdepth[name],
                                                                      sample[analysistype].standarddev[name])
                        # Add a newline after each sample
                        data += '\n'
                    # Add a newline if the sample did not have any gene hits
                    else:
                        data += '\n'
            # Write the header and data to file
            report.write(header)
            report.write(data)

    def genusspecific(self, analysistype='genesippr'):
        """
        Creates simplified genus-specific reports. Instead of the % ID and the fold coverage, a simple +/- scheme is
        used for presence/absence
        :param analysistype: The variable to use when accessing attributes in the metadata object
        """
        # Dictionary containing genera of interest, and the probes in the database
        genedict = {'Escherichia': ['eae', 'O26', 'O45', 'O103', 'O111', "O121", 'O145', 'O157', 'VT1', 'VT2',
                                    'VT2f', 'uidA'],
                    'Listeria': ['hlyA', 'IGS', 'inlJ'],
                    'Salmonella': ['invA', 'stn']}
        # Dictionary to store all the output strings
        results = dict()
        for genus, genelist in genedict.items():
            # Initialise the dictionary with the appropriate genus
            results[genus] = str()
            for sample in self.runmetadata.samples:
                try:
                    # Find the samples that match the current genus - note that samples with multiple hits will be
                    # represented in multiple outputs
                    if genus in sample[analysistype].targetgenera:
                        # Populate the results string with the sample name
                        results[genus] += '{},'.format(sample.name)
                        # Iterate through all the genes associated with this genus. If the gene is in the current
                        # sample, add a + to the string, otherwise, add a -
                        for gene in genelist:
                            if gene.lower() in [target[0].lower() for target in sample[analysistype].results.items()]:
                                results[genus] += '+,'
                            # Special case for eae: if either eae or eae_1 are present, add a +
                            elif gene.lower() + '_1' in [target[0].lower() for target in
                                                         sample[analysistype].results.items()]:
                                results[genus] += '+,'
                            else:
                                results[genus] += '-,'
                        results[genus] += '\n'
                # If the sample is missing the targetgenera attribute, then it is ignored for these reports
                except KeyError:
                    pass
        # Create and populate the genus-specific reports
        for genus, resultstring in results.items():
            # Only create the report if there are results for the current genus
            if resultstring:
                with open(os.path.join(self.reportpath, '{}_genesippr.csv'.format(genus)), 'w') as genusreport:
                    # Write the header to the report - Strain plus add the genes associated with the genus
                    genusreport.write('Strain,{}\n'.format(','.join(genedict[genus])))
                    # Write the results to the report
                    genusreport.write(resultstring)

    def gdcsreporter(self, analysistype='GDCS'):
        """
        Creates a report of the GDCS results
        :param analysistype: The variable to use when accessing attributes in the metadata object
        """
        printtime('Creating {} report'.format(analysistype), self.starttime, output=self.portallog)
        # Initialise list to store all the GDCS genes, and genera in the analysis
        gdcs = list()
        genera = list()
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                sample[analysistype].createreport = True
                # Determine which genera are present in the analysis
                if sample.general.closestrefseqgenus not in genera:
                    genera.append(sample.general.closestrefseqgenus)
                try:
                    # Add all the GDCS genes to the list
                    for gene in sorted(sample[analysistype].faidict):
                        if gene not in gdcs:
                            gdcs.append(gene)
                except KeyError:
                    sample[analysistype].createreport = False
            else:
                sample[analysistype].createreport = False
        header = 'Strain,Genus,Matches,MeanCoverage,Pass/Fail,{},\n'.format(','.join(gdcs))
        data = str()
        with open(os.path.join(self.reportpath, '{}.csv'.format(analysistype)), 'w') as report:
            # Sort the samples in the report based on the closest refseq genus e.g. all samples with the same genus
            # will be grouped together in the report
            for genus in genera:
                for sample in self.runmetadata.samples:
                    if sample.general.closestrefseqgenus == genus:
                        if sample[analysistype].createreport:
                            sample[analysistype].totaldepth = list()
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
                                if gene not in [result for result in sample[analysistype].faidict]:
                                    specific += 'X,'
                                else:
                                    try:
                                        # Report the necessary information for each gene result
                                        identity = sample[analysistype].results[gene]
                                        specific += '{}% ({} +/- {}),'\
                                            .format(identity, sample[analysistype].avgdepth[gene],
                                                    sample[analysistype].standarddev[gene])
                                        sample[analysistype].totaldepth.append(
                                            float(sample[analysistype].avgdepth[gene]))
                                        count += 1
                                    # If the gene was missing from the results attribute, add a - to the cell
                                    except KeyError:
                                        sample.general.incomplete = True
                                        specific += '-,'
                            # Calculate the mean depth of the genes and the standard deviation
                            sample[analysistype].mean = numpy.mean(sample[analysistype].totaldepth)
                            sample[analysistype].stddev = numpy.std(sample[analysistype].totaldepth)
                            # Determine whether the sample pass the necessary quality criteria:
                            # Pass, all GDCS, mean coverage greater than 20X coverage;
                            # ?: Indeterminate value;
                            # -: Fail value
                            if count == len(sample[analysistype].faidict):
                                if sample[analysistype].mean > 20:
                                    quality = '+'
                                else:
                                    quality = '?'
                                    sample.general.incomplete = True
                            else:
                                quality = '-'
                                sample.general.incomplete = True
                            # Add the count, mean depth with standard deviation, the pass/fail determination,
                            #  and the total number of GDCS genes as well as the results
                            data += '{hits}/{total},{mean} +/- {std},{fail},{gdcs}\n'\
                                .format(hits=str(count),
                                        total=len(sample[analysistype].faidict),
                                        mean='{:.2f}'.format(sample[analysistype].mean),
                                        std='{:.2f}'.format(sample[analysistype].stddev),
                                        fail=quality,
                                        gdcs=specific)
                        # Any samples with a best assembly of 'NA' are considered incomplete.
                        else:
                            data += '{},{},,,-\n'.format(sample.name, sample.general.closestrefseqgenus)
                            sample.general.incomplete = True
                    elif sample.general.closestrefseqgenus == 'NA':
                        data += '{}\n'.format(sample.name)
                        sample.general.incomplete = True
            # Write the header and data to file
            report.write(header)
            report.write(data)

    def sixteensreporter(self, analysistype='sixteens_full'):
        """
        Creates a report of the results
        :param analysistype: The variable to use when accessing attributes in the metadata object
        """
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # Initialise the header and data strings
        header = 'Strain,Gene,PercentIdentity,Genus,FoldCoverage\n'
        data = ''
        with open(os.path.join(self.reportpath, analysistype + '.csv'), 'w') as report:
            with open(os.path.join(self.reportpath, analysistype + '_sequences.fa'), 'w') as sequences:
                for sample in self.runmetadata.samples:
                    try:
                        # Select the best hit of all the full-length 16S genes mapped
                        sample[analysistype].besthit = sorted(sample[analysistype].results.items(),
                                                              key=operator.itemgetter(1), reverse=True)[0][0]
                        # Add the sample name to the data string
                        data += sample.name + ','
                        # Find the record that matches the best hit, and extract the necessary values to be place in the
                        # data string
                        for name, identity in sample[analysistype].results.items():
                            if name == sample[analysistype].besthit:
                                data += '{},{},{},{}\n'.format(name, identity, sample[analysistype].genus,
                                                               sample[analysistype].avgdepth[name])
                                # Create a FASTA-formatted sequence output of the 16S sequence
                                record = SeqRecord(Seq(sample[analysistype].sequences[name],
                                                       IUPAC.unambiguous_dna),
                                                   id='{}_{}'.format(sample.name, '16S'),
                                                   description='')
                                SeqIO.write(record, sequences, 'fasta')
                    except (KeyError, IndexError):
                        data += '{}\n'.format(sample.name)
            # Write the results to the report
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
        # Create the reports
        self.reporter()
        self.genusspecific()
        self.sixteensreporter()
        self.gdcsreporter()

    def __init__(self, inputobject):
        self.starttime = inputobject.starttime
        self.samples = inputobject.samples
        self.completemetadata = inputobject.completemetadata
        self.path = inputobject.path
        self.analysescomplete = inputobject.analysescomplete
        self.reportpath = inputobject.reportpath
        self.runmetadata = inputobject.runmetadata
        try:
            self.portallog = inputobject.portallog
        except AttributeError:
            self.portallog = ''

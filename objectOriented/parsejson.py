#!/usr/bin/env python
import json

from SPAdesPipeline.OLCspades.accessoryFunctions import *


class ParseJson(object):
    def findjson(self):
        """
        Find the json files and create objects to store metadata
        """
        from glob import glob
        printtime('Finding JSON reports and creating python data objects', self.start)
        # Get a list of the json files
        jsonfiles = glob('{}*/*.json'.format(self.path))
        for jsonfile in sorted(jsonfiles):
            # Create the object to store the metadata information
            metadata = MetadataObject()
            metadata.general = GenObject()
            # Populate the attributes of the object as necessary
            metadata.general.jsonfile = jsonfile
            metadata.name = os.path.basename(jsonfile).split('_')[0]
            metadata.general.outdir = os.path.split(jsonfile)[0]
            self.samples.append(metadata)

    def parse(self):
        """
        Parse the json file and create a fasta file of all the new, non-redundant, alleles
        """
        # Import BioPython modules for creating a nice output
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Set the name and path of the file containing the new alleles
        newalleles = '{}newalleles.fa'.format(self.path)
        printtime('Parsing results and creating file of new alleles. Path of file: {}'.format(newalleles), self.start)
        # Initialise lists
        alleles = list()
        alleletuple = list()
        for sample in self.samples:
            # Open the json file
            with open(sample.general.jsonfile, 'rb') as jsondata:
                # Load the data from the json file into a dictionary
                jsondict = json.load(jsondata)
                # In this script, I'm only looking for new alleles within the 'rmlst' analyses
                for gene, allele in jsondict['rmlst']['newalleles'].items():
                    # Only add alleles if they are not already in the list
                    if allele not in alleles:
                        alleles.append(allele)
                        # Create a tuple of the gene name and the allele sequence
                        alleletuple.append((gene, allele))
        # Open the file to write the new alleles
        with open(newalleles, 'wb') as allelefile:
            # Iterate through the allele data sorted on gene name
            for alleledata in sorted(alleletuple):
                # Create a sequence record using BioPython
                fasta = SeqRecord(Seq(alleledata[1]),
                                  # Without this, the header will be improperly formatted
                                  description='',
                                  # Use >:definitionline as the header
                                  id=alleledata[0])
                # Use the SeqIO module to properly format the new sequence record
                SeqIO.write(fasta, allelefile, "fasta")

    def __init__(self, args, starttime):
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.start = starttime
        self.samples = list()
        # Call the necessary modules
        self.findjson()
        self.parse()


if __name__ == '__main__':
    import time
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser

    # Parser for arguments
    parser = ArgumentParser(description='Find new alleles written to JSON metadata files')
    parser.add_argument('path',
                        help='Specify input directory. Note that this path must be the same as '
                             'the sequence path specified in geneSippr')

    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    start = time.time()

    # Run the script
    ParseJson(arguments, start)

    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m'

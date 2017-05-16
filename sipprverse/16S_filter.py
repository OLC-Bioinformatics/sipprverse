#!/usr/bin/env python
from sipprcommon.accessoryfunctions.accessoryFunctions import *
__author__ = 'adamkoziol'


class Filter(object):

    def filter(self):
        """Filter all desired sequences into a .fasta file"""
        from Bio import SeqIO
        # A list to store all the records from the organism(s) of interest
        targets = list()
        # Set the name of the filtered assembly file
        filteredfile = os.path.join(self.output, '{}_filtered.fasta'
                                    .format(os.path.basename(self.targetfile).split('.')[0]))
        # Find all the records corresponding to the desired genus, sequence binomials
        for genus, species in self.organisms.items():
            printtime('Processing {} {}'.format(genus, species), self.starttime)
            for record in SeqIO.parse(open(self.targetfile, "rU"), "fasta"):
                # Include only sequences that contain both the genus and the species
                if genus in record.description and species in record.description:
                    # Add this record to our list
                    print(record.description)
                    targets.append(record)
        # Only create the file if there are contigs over 1000 bp
        if targets:
            # Open the filtered assembly file
            with open(filteredfile, 'wb') as formatted:
                # Write the records in the list to the file
                SeqIO.write(targets, formatted, 'fasta')

    def customfilter(self):
        """Filter all desired sequences into a .fasta file"""
        from Bio import SeqIO
        # A list to store all the records of interest
        targets = list()
        # Set the name of the filtered assembly file
        filteredfile = os.path.join(self.output, '{}_{}_filtered.fasta'
                                    .format(os.path.basename(self.targetfile).split('.')[0], '_'.join(self.custom)))
        # Find all the records with the custom string in the header
        printtime('Searching for sequences with {} in the header'.format(self.custom), self.starttime)
        for record in SeqIO.parse(open(self.targetfile, "rU"), "fasta"):
            # Include only sequences that contain the custom string, and do not contain any of the strings to ignore
            if all(custom.lower() in record.description.lower() for custom in self.custom) \
                    and not any(ignore.lower() in record.description.lower() for ignore in self.ignore):
                # Add this record to our list
                print(record.description)
                targets.append(record)
        # Only create the file if there are contigs over 1000 bp
        if targets:
            # Open the filtered assembly file
            with open(filteredfile, 'wb') as formatted:
                # Write the records in the list to the file
                SeqIO.write(targets, formatted, 'fasta')

    def __init__(self, args):
        """
        :param args: command line arguments
        """
        # Initialise variables
        self.starttime = args.start
        # Define variables based on supplied arguments
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.targetfile = os.path.join(self.path, args.targetfile)
        #
        assert os.path.isfile(self.targetfile), u'Target file does not exist in the path {0!r:s}' \
            .format(self.targetfile)
        self.output = os.path.join(self.path, 'filtered')
        make_path(self.output)
        try:
            self.custom = args.custom.split(',')
        except AttributeError:
            self.custom = ''
        if not self.custom:
            # Create a dict of genus: species of the organisms to use by creating a list of the genus_species entries
            # with a list comprehension splitting on commas. Split the genus from the species by mapping the split on
            # underscores to this list of genus_species. Convert the list of lists created by the map function to a
            # dictionary
            self.organisms = dict(map(lambda genus: genus.split('_'), [org for org in args.organisms.split(',')]))
            # Filter the input file with the organism dictionary
            self.filter()
        # Run the custom string searching/ignoring method
        else:
            try:
                self.ignore = args.ignore.split(',')
            # If no ignore string is provided, use an empty string
            except AttributeError:
                self.ignore = ''
            self.customfilter()


if __name__ == '__main__':
    import time
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Extract 16S sequences from a .fasta file for desired organisms. '
                                        'Alternatively, extract any sequence from a .fasta file that contains a '
                                        'custom string in the header, and, optionally, does not contain a different '
                                        'string in the header')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-t', '--targetfile',
                        required=True,
                        help='Name of sixteenS target file to process. Note that this file must be within the supplied '
                             'path, and that this file must have the target organism name in the header. Or, if using '
                             '--custom, it can have any FASTA format')
    parser.add_argument('-o', '--organisms',
                        default='Escherichia_coli,Salmonella_enterica,Listeria_monocytogenes',
                        help='Supply a comma-separated list of organisms to use. The default is '
                             '"Escherichia_coli,Salmonella_enterica,Listeria_monocytogenes". Note that the genus and '
                             'species of the organisms are separated by an underscore')
    parser.add_argument('-c', '--custom',
                        help='Comma-separated list. Extract sequences containing a custom string e.g.stx1,:a will '
                             'extract sequences containing both strings (stx1B:3:AB083044:a)')
    parser.add_argument('-i', '--ignore',
                        help='Comma-separated list. Used with the --custom flag. Will exclude sequences that '
                             'contain the string e.g. :f,:g will cause anything that is matched to be ignored if it '
                             'has :f or :g in the header')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.start = time.time()

    # Run the script
    Filter(arguments)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.start) + '\033[0m')

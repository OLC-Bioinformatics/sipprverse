#!/usr/bin/env python
import subprocess
import time
import os
from MLSTsippr.sipprmlst import MLSTmap
from sipprCommon.objectprep import Objectprep
from accessoryFunctions.accessoryFunctions import printtime, make_dict, dotter, make_path
from accessoryFunctions.metadataprinter import MetadataPrinter
from collections import defaultdict
import operator
__author__ = 'adamkoziol'


class GeneSippr(object):

    def runner(self):
        """
        Run the necessary methods in the correct order
        """
        printtime('Starting {} analysis pipeline'.format(self.analysistype), self.starttime)
        # Create the objects to be used in the analyses (if required)
        general = None
        for sample in self.runmetadata.samples:
            general = getattr(sample, 'general')
        if general is None:
            # Create the objects to be used in the analyses
            objects = Objectprep(self)
            objects.objectprep()
            self.runmetadata = objects.samples
        # Run the analyses
        MLSTmap(self, self.analysistype, self.cutoff)
        # Create the reports
        self.reporter()
        for sample in self.runmetadata.samples:
            # Remove large attributes from the object
            try:
                delattr(sample[self.analysistype], 'profiledata')
                delattr(sample[self.analysistype], 'allelenames')
                delattr(sample[self.analysistype], 'alleles')
                delattr(sample[self.analysistype], 'faidict')
                delattr(sample[self.analysistype], 'gaplocations')
                delattr(sample[self.analysistype], 'maxcoverage')
                delattr(sample[self.analysistype], 'mincoverage')
                delattr(sample[self.analysistype], 'resultsgap')
                delattr(sample[self.analysistype], 'resultssnp')
                delattr(sample[self.analysistype], 'snplocations')
                delattr(sample[self.analysistype], 'standarddev')
                delattr(sample[self.analysistype], 'avgdepth')
            except KeyError:
                pass
        # Print the metadata to a .json file
        MetadataPrinter(self)

    def reporter(self):
        """
        Runs the necessary methods to parse raw read outputs
        """
        printtime('Preparing reports', self.starttime)
        # Populate self.plusdict in order to reuse parsing code from an assembly-based method
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                for gene in sample[self.analysistype].allelenames:
                    for allele, percentidentity in sample[self.analysistype].results.items():
                        if gene in allele:
                            # Split the allele number from the gene name using the appropriate delimiter
                            if '_' in allele:
                                splitter = '_'
                            elif '-' in allele:
                                splitter = '-'
                            else:
                                splitter = ''
                            # Create the plusdict dictionary as in the assembly-based (r)MLST method. Allows all the
                            # parsing and sequence typing code to be reused.
                            try:
                                self.plusdict[sample.name][gene][allele.split(splitter)[1]][percentidentity] \
                                    = sample[self.analysistype].avgdepth[allele]
                            except IndexError:
                                pass
        self.profiler()
        self.sequencetyper()
        self.mlstreporter()

    def profiler(self):
        """Creates a dictionary from the profile scheme(s)"""
        printtime('Loading profiles', self.starttime)
        from csv import DictReader
        # Initialise variables
        profiledata = defaultdict(make_dict)
        profileset = set()
        genedict = dict()
        # Find all the unique profiles to use with a set
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                if sample[self.analysistype].profile != 'NA':
                    profileset.add(sample[self.analysistype].profile)

        # Extract the profiles for each set
        for sequenceprofile in profileset:
            # Clear the list of genes
            genelist = list()
            for sample in self.runmetadata.samples:
                if sample.general.bestassemblyfile != 'NA':
                    if sequenceprofile == sample[self.analysistype].profile:
                        genelist = [allele for allele in sample[self.analysistype].alleles]
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
                        try:
                            profiledata[sequenceprofile][row['rST']][gene] = row[gene]
                        except KeyError:
                            raise
            # Add the gene list to a dictionary
            genedict[sequenceprofile] = sorted(genelist)
            # Add the profile data, and gene list to each sample
            for sample in self.runmetadata.samples:
                if sample.general.bestassemblyfile != 'NA':
                    if sequenceprofile == sample[self.analysistype].profile:
                        # Populate the metadata with the profile data
                        sample[self.analysistype].profiledata = profiledata[sample[self.analysistype].profile]
                        dotter()

    def sequencetyper(self):
        """Determines the sequence type of each strain based on comparisons to sequence type profiles"""
        printtime('Performing sequence typing', self.starttime)
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                if type(sample[self.analysistype].allelenames) == list:
                    # Initialise variables
                    header = 0
                    # Iterate through the genomes
                    # for sample in self.runmetadata.samples:
                    genome = sample.name
                    # Initialise self.bestmatch[genome] with an int that will eventually be replaced by the # of matches
                    self.bestmatch[genome] = defaultdict(int)
                    if sample[self.analysistype].profile != 'NA':
                        # Create the profiledata variable to avoid writing self.profiledata[self.analysistype]
                        profiledata = sample[self.analysistype].profiledata
                        # For each gene in plusdict[genome]
                        for gene in sample[self.analysistype].allelenames:
                            # Clear the appropriate count and lists
                            multiallele = list()
                            multipercent = list()
                            # Go through the alleles in plusdict
                            for allele in self.plusdict[genome][gene]:
                                percentid = list(self.plusdict[genome][gene][allele].keys())[0]
                                # "N" alleles screw up the allele splitter function
                                if allele != "N":
                                    # Use the alleleSplitter function to get the allele number
                                    # allelenumber, alleleprenumber = allelesplitter(allele)
                                    # Append as appropriate - alleleNumber is treated as an integer for proper sorting
                                    multiallele.append(int(allele))
                                    multipercent.append(percentid)
                                # If the allele is "N"
                                else:
                                    # Append "N" and a percent identity of 0
                                    multiallele.append("N")
                                    multipercent.append(0)
                            # Populate self.bestdict with genome, gene, alleles joined with a space (this was made like
                            # this because allele is a list generated by the .iteritems() above
                            try:
                                self.bestdict[genome][gene][" ".join(str(allele)
                                                                     for allele in sorted(multiallele))] = \
                                    multipercent[0]
                            except IndexError:
                                self.bestdict[genome][gene]['NA'] = 0
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
                                for allele, percentid in self.bestdict[genome][gene].items():
                                    # If the allele in the query genome matches the allele in the reference profile, add
                                    # the result to the bestmatch dictionary. Genes with multiple alleles were sorted
                                    # the same, strings with multiple alleles will match: 10 692 will never be 692 10
                                    if allele == sortedrefallele and float(percentid) == 100.00:
                                        # Increment the number of matches to each profile
                                        self.bestmatch[genome][sequencetype] += 1
                                    # Special handling of BACT000060 and BACT000065 genes for E. coli and BACT000014
                                    # for Listeria. When the reference profile has an allele of 'N', and the query
                                    # allele doesn't, set the allele to 'N', and count it as a match
                                    elif gene == 'BACT000060' or gene == 'BACT000065' or gene == 'BACT000014':
                                        if sortedrefallele == 'N' and allele != 'N':
                                            # Increment the number of matches to each profile
                                            self.bestmatch[genome][sequencetype] += 1
                                    elif allele == sortedrefallele and sortedrefallele == 'N':
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
                            for sequencetype, matches in self.bestmatch[genome].items():
                                if matches == sortedmatches:
                                    for gene in profiledata[sequencetype]:
                                        # Populate resultProfile with the genome, best match to profile, # of matches
                                        # to the profile, gene, query allele(s), reference allele(s), and % identity
                                        self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                            list(self.bestdict[genome][gene]
                                                 .keys())[0]] = str(list(self.bestdict[genome][gene].values())[0])
                                    sample[self.analysistype].sequencetype = sequencetype
                                    sample[self.analysistype].matchestosequencetype = matches
                        # If there are fewer matches than the total number of genes in the typing scheme
                        elif 0 < int(sortedmatches) < header:
                            mismatches = []
                            # Iterate through the sequence types and the number of matches in bestDict for each genome
                            for sequencetype, matches in self.bestmatch[genome].items():
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
                                        if self.analysistype == 'mlst':
                                            self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                                list(self.bestdict[genome][gene]
                                                     .keys())[0]] = str(list(self.bestdict[genome][gene].values())[0])
                                            sample[self.analysistype].mismatchestosequencetype = mismatches
                                            sample[self.analysistype].sequencetype = sequencetype
                                            sample[self.analysistype].matchestosequencetype = matches
                                        else:
                                            self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                                list(self.bestdict[genome][gene].keys())[0]] \
                                                = str(list(self.bestdict[genome][gene].values())[0])
                                            #
                                            if sortedrefallele != list(self.bestdict[sample.name][gene].keys())[0]:
                                                mismatches.append(
                                                    ({gene: ('{} ({})'.format(list(self.bestdict[sample.name][gene]
                                                                                   .keys())[0], sortedrefallele))}))
                        elif sortedmatches == 0:
                            for gene in sample[self.analysistype].allelenames:
                                # Populate the results profile with negative values for sequence type and sorted matches
                                self.resultprofile[genome]['NA'][sortedmatches][gene]['NA'] = 0
                            # Add the new profile to the profile file (if the option is enabled)
                            sample[self.analysistype].sequencetype = 'NA'
                            sample[self.analysistype].matchestosequencetype = 'NA'
                            sample[self.analysistype].mismatchestosequencetype = 'NA'
                        else:
                            sample[self.analysistype].matchestosequencetype = 'NA'
                            sample[self.analysistype].mismatchestosequencetype = 'NA'
                            sample[self.analysistype].sequencetype = 'NA'
                        dotter()
                else:
                    sample[self.analysistype].matchestosequencetype = 'NA'
                    sample[self.analysistype].mismatchestosequencetype = 'NA'
                    sample[self.analysistype].sequencetype = 'NA'

            else:
                sample[self.analysistype].matchestosequencetype = 'NA'
                sample[self.analysistype].mismatchestosequencetype = 'NA'
                sample[self.analysistype].sequencetype = 'NA'

    def mlstreporter(self):
        """ Parse the results into a report"""
        printtime('Writing reports', self.starttime)
        # Initialise variables
        combinedrow = str()
        reportdirset = set()
        # Populate a set of all the report directories to use. A standard analysis will only have a single report
        # directory, while pipeline analyses will have as many report directories as there are assembled samples
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                # Ignore samples that lack a populated reportdir attribute
                if sample[self.analysistype].reportdir != 'NA':
                    make_path(sample[self.analysistype].reportdir)
                    # Add to the set - I probably could have used a counter here, but I decided against it
                    reportdirset.add(sample[self.analysistype].reportdir)
        # Create a report for each sample from :self.resultprofile
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                if sample[self.analysistype].reportdir != 'NA':
                    if type(sample[self.analysistype].allelenames) == list:
                        # Populate the header with the appropriate data, including all the genes in the list of targets
                        row = 'Strain,Genus,SequenceType,Matches,{},\n' \
                            .format(','.join(sorted(sample[self.analysistype].allelenames)))
                        # Set the seq counter to 0. This will be used when a sample has multiple best sequence types.
                        # The sample name will not be written on subsequent rows in order to make the report clearer
                        seqcount = 0
                        # Iterate through the best sequence types for the sample
                        for seqtype in self.resultprofile[sample.name]:
                            sample[self.analysistype].sequencetype = seqtype
                            # The number of matches to the profile
                            sample[self.analysistype].matches = list(self.resultprofile[sample.name][seqtype].keys())[0]
                            # If this is the first of one or more sequence types, include the sample name
                            if seqcount == 0:
                                row += '{},{},{},{},'.format(sample.name, sample.general.referencegenus, seqtype,
                                                             sample[self.analysistype].matches)
                            # Otherwise, skip the sample name
                            else:
                                row += ',,{},{},'.format(seqtype, sample[self.analysistype].matches)
                            # Iterate through all the genes present in the analyses for the sample
                            for gene in sorted(sample[self.analysistype].allelenames):
                                refallele = sample[self.analysistype].profiledata[seqtype][gene]
                                # Set the allele and percent id from the dictionary's keys and values, respectively
                                allele = \
                                    list(self.resultprofile[sample.name][seqtype][sample[self.analysistype].matches]
                                         [gene].keys())[0]
                                percentid = \
                                    list(self.resultprofile[sample.name][seqtype][sample[self.analysistype].matches]
                                         [gene].values())[0]
                                try:
                                    if refallele and refallele != allele:
                                        if 0 < float(percentid) < 100:
                                            row += '{} ({:.2f}%),'.format(allele, float(percentid))
                                        else:
                                            row += '{} ({}),'.format(allele, refallele)
                                    else:
                                        # Add the allele and % id to the row (only add the % identity if it is not 100%)
                                        if 0 < float(percentid) < 100:
                                            row += '{} ({:.2f}%),'.format(allele, float(percentid))
                                        else:
                                            row += '{},'.format(allele)
                                    self.referenceprofile[sample.name][gene] = allele
                                except ValueError:
                                    pass
                            # Add a newline
                            row += '\n'
                            # Increment the number of sequence types observed for the sample
                            seqcount += 1
                        combinedrow += row
                        # If the length of the # of report directories is greater than 1 (script is being run as part of
                        # the assembly pipeline) make a report for each sample
                        if self.pipeline:
                            # Open the report
                            with open(os.path.join(sample[self.analysistype].reportdir,
                                                   '{}_{}.csv'.format(sample.name, self.analysistype)), 'w') as report:
                                # Write the row to the report
                                report.write(row)
                dotter()
            # Create the report folder
            make_path(self.reportpath)
            # Create the report containing all the data from all samples
            if self.pipeline:
                with open(os.path.join(self.reportpath,  '{}.csv'.format(self.analysistype)), 'w') \
                        as combinedreport:
                    # Write the results to this report
                    combinedreport.write(combinedrow)
            else:
                with open(os.path.join(self.reportpath, '{}_{:}.csv'.format(
                        self.analysistype, time.strftime("%Y.%m.%d.%H.%M.%S"))), 'w') as combinedreport:
                    # Write the results to this report
                    combinedreport.write(combinedrow)

    def __init__(self, args, pipelinecommit, startingtime, scriptpath, analysistype, cutoff, pipeline):
        """
        :param args: command line arguments
        :param pipelinecommit: pipeline commit or version
        :param startingtime: time the script was started
        :param scriptpath: home path of the script
        :param analysistype: name of the analysis being performed - allows the program to find databases
        :param cutoff: percent identity cutoff for matches
        :param pipeline: boolean of whether this script needs to run as part of a particular assembly pipeline
        """
        import multiprocessing
        # Initialise variables
        self.commit = str(pipelinecommit)
        self.starttime = startingtime
        self.homepath = scriptpath
        # Define variables based on supplied arguments
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        try:
            self.sequencepath = os.path.join(args.sequencepath, '')
        except AttributeError:
            self.sequencepath = self.path
        assert os.path.isdir(self.sequencepath), u'Sequence path  is not a valid directory {0!r:s}' \
            .format(self.sequencepath)
        try:
            self.targetpath = os.path.join(args.reffilepath, analysistype)
        except AttributeError:
            self.targetpath = os.path.join(args.targetpath, '')
        assert os.path.isdir(self.targetpath), u'Target path is not a valid directory {0!r:s}' \
            .format(self.targetpath)
        self.reportpath = os.path.join(self.path, 'reports')
        try:
            self.bcltofastq = args.bcltofastq
        except AttributeError:
            self.bcltofastq = False
        try:
            self.miseqpath = args.miseqpath
        except AttributeError:
            self.miseqpath = str()
        try:
            self.miseqfolder = args.miseqfolder
        except AttributeError:
            self.miseqfolder = str()
        try:
            self.fastqdestination = args.fastqdestination
        except AttributeError:
            self.fastqdestination = str()
        try:
            self.forwardlength = args.forwardlength
        except AttributeError:
            self.forwardlength = 'full'
        try:
            self.reverselength = args.reverselength
        except AttributeError:
            self.reverselength = 'full'
        self.numreads = 2 if self.reverselength != 0 else 1
        self.customsamplesheet = args.customsamplesheet
        # Set the custom cutoff value
        self.cutoff = float(cutoff)
        self.logfile = args.logfile
        try:
            self.averagedepth = int(args.averagedepth)
        except AttributeError:
            self.averagedepth = 10
        try:
            self.copy = args.copy
        except AttributeError:
            self.copy = False
        self.runmetadata = args.runmetadata
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        try:
            self.cpus = int(args.cpus)
        except AttributeError:
            self.cpus = multiprocessing.cpu_count()
        try:
            self.threads = int(self.cpus / len(self.runmetadata.samples)) if self.cpus / len(self.runmetadata.samples) \
                                                                             > 1 else 1
        except TypeError:
            self.threads = self.cpus
        self.taxonomy = {'Escherichia': 'coli', 'Listeria': 'monocytogenes', 'Salmonella': 'enterica'}
        #
        self.pipeline = pipeline
        if analysistype.lower() == 'mlst':
            self.analysistype = 'mlst'
        elif analysistype.lower() == 'rmlst':
            self.analysistype = 'rmlst'
        else:
            import sys
            sys.stderr.write('Please ensure that you specified a valid option for the analysis type. You entered {}. '
                             'The only acceptable options currently are mlst and rmlst.'.format(args.analysistype))
            quit()
        self.plusdict = defaultdict(make_dict)
        self.bestdict = defaultdict(make_dict)
        self.bestmatch = defaultdict(int)
        self.mlstseqtype = defaultdict(make_dict)
        self.resultprofile = defaultdict(make_dict)
        self.referenceprofile = defaultdict(make_dict)
        # Run the analyses
        self.runner()


if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Parser for arguments
    parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    parser.add_argument('-t', '--targetpath',
                        required=True,
                        help='Path of target files to process.')
    parser.add_argument('-n', '--numthreads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-b', '--bcl2fastq',
                        action='store_true',
                        help='Optionally run bcl2fastq on an in-progress Illumina MiSeq run. Must include:'
                             'miseqpath, and miseqfolder arguments, and optionally readlengthforward, '
                             'readlengthreverse, and projectName arguments.')
    parser.add_argument('-m', '--miseqpath',
                        help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', '--miseqfolder',
                        help='Name of the folder containing MiSeq run data')
    parser.add_argument('-d', '--fastqdestination',
                        help='Optional folder path to store .fastq files created using the fastqCreation module. '
                             'Defaults to path/miseqfolder')
    parser.add_argument('-r1', '--forwardlength',
                        default='full',
                        help='Length of forward reads to use. Can specify "full" to take the full length of '
                             'forward reads specified on the SampleSheet')
    parser.add_argument('-r2', '--reverselength',
                        default='full',
                        help='Length of reverse reads to use. Can specify "full" to take the full length of '
                             'reverse reads specified on the SampleSheet')
    parser.add_argument('-c', '--customsamplesheet',
                        help='Path of folder containing a custom sample sheet (still must be named "SampleSheet.csv")')
    parser.add_argument('-P', '--projectName',
                        help='A name for the analyses. If nothing is provided, then the "Sample_Project" field '
                             'in the provided sample sheet will be used. Please note that bcl2fastq creates '
                             'subfolders using the project name, so if multiple names are provided, the results '
                             'will be split as into multiple projects')
    parser.add_argument('-D', '--detailedReports',
                        action='store_true',
                        help='Provide detailed reports with percent identity and depth of coverage values '
                             'rather than just "+" for positive results')
    parser.add_argument('-u', '--customcutoffs',
                        default=1.0,
                        help='Custom cutoff values')
    parser.add_argument('-a', '--analysistype',
                        required=True,
                        help='Specify analysis type: mlst or rmlst')
    parser.add_argument('-C', '--copy',
                        action='store_true',
                        help='Normally, the program will create symbolic links of the files into the sequence path, '
                             'however, the are occasions when it is necessary to copy the files instead')

    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.pipeline = False
    arguments.logfile = os.path.join(arguments.path, 'logfile')
    # Define the start time
    start = time.time()

    # Run the script
    GeneSippr(arguments, commit, start, homepath, arguments.analysistype, arguments.customcutoffs, arguments.pipeline)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')

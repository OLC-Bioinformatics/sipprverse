#!/usr/bin/env python 3
from sipprcommon.accessoryfunctions.accessoryFunctions import *
from Bio import SeqIO
from ete3 import NCBITaxa
import shutil
import pickle
__author__ = 'adamkoziol'


class Filter(object):

    def runner(self):
        """

        """
        #
        self.databasereduce()
        self.databaseparse()
        #
        # self.paths()
        #
        self.cluster()
        self.clustertargets()
        # self.align()
        #
        # self.probefinder()
        #
        # self.probes()
        #
        # self.targets()

    def databasereduce(self):
        """
        Create a reduced taxonomy file containing only the necessary taxonomy IDs from the nucl_gb.accession2taxid
        that are present in the file. Helpful if you are planning on parsing the database more than once
        """
        if not os.path.isfile(self.pathobject):
            # List to store all the accessions of the files in your database
            records = set()
            printtime('Reading in {} records'.format(self.analysistype), self.starttime)
            # Index the multifasta file, and extract all the accessions
            self.recorddict = SeqIO.index(self.targetfile, 'fasta')
            if self.analysistype == 'sixteens':
                for record in self.recorddict:
                    accession = record.split('|')[-2]
                    self.linkingdict[accession] = record
                    records.add(accession)
            else:
                with open(self.targetfile, 'rU') as targetfile:
                    for record in SeqIO.parse(targetfile, 'fasta'):
                        accession = record.description.split()[1]
                        self.linkingdict[accession] = record.id
                        records.add(accession)
            if not os.path.isfile(self.reducedtaxidfile):
                printtime('Loading taxonomy information', self.starttime)
                # Dictionary to store all the lines that match an accession number from :records
                reduced = dict()
                # Use nucl_gb.accession2taxid to determine the taxonomy ID for each accession
                with open(self.taxidfile, 'r') as taxidfile:
                    for line in taxidfile:
                        # The taxid file has the following format:
                        '''
                        accession	accession.version	taxid	gi
                        A00002	A00002.1	9913	2
                        '''
                        if self.analysistype == 'sixteens':
                            # The accession.version is what needs to match the accessions in :records
                            _, accession, taxid, _ = line.split('\t')
                        else:
                            accession, _, taxid, _ = line.split('\t')
                        # If the accession.version is in the list, add the accession.version: taxid to the dictionary
                        if accession in records:
                            reduced[accession] = taxid
                printtime('Creating reduced taxonomy file', self.starttime)
                # Use pickle to serialise the dictionary of filtered taxids
                with open(self.reducedtaxidfile, 'wb') as reducedfile:
                    pickle.dump(reduced, reducedfile)

    def databaseparse(self):
        """
        Parses filtered taxonomy IDs, and sorts the files based on taxonomic level
        """
        printtime('Splitting {} records into taxonomic groups'.format(self.analysistype), self.starttime)
        if not os.path.isfile(self.pathobject):
            # Create an NCBITaxa object to be used to query the taxonomy database
            ncbi = NCBITaxa()
            # Remove the output directory - as all the output files are always appended, running the script multiple
            # times will create large, useless files
            try:
                shutil.rmtree(self.output)
                shutil.rmtree(self.genuspath)
            # If the folder doesn't exist, ignore
            except FileNotFoundError:
                pass
            make_path(self.output)
            bacteriafile = os.path.join(self.path, 'bacterial_filtered.fa')

            make_path(self.genuspath)
            try:
                os.remove(bacteriafile)
            except OSError:
                pass
            # Load the pickled filtered taxids into a dictionary
            with open(self.reducedtaxidfile, 'rb') as taxidfile:
                taxiddict = pickle.load(taxidfile)
                for accession, taxid in taxiddict.items():
                    # Allow for missing tax IDs
                    try:
                        # Use NCBITaxa to get the lineage data for each taxid
                        lineage = ncbi.get_lineage(taxid)
                        # Convert the lineage taxids to a dictionary of taxonomy names
                        names = ncbi.get_taxid_translator(lineage)
                        # List comprehension to convert dictionary to list - to allow sorting
                        taxonomylist = [names[tax] for tax in lineage]
                        # Only process bacteria
                        if 'Bacteria' in taxonomylist:
                            # Reset the path for each accession
                            path = self.output
                            # As the path is being updated below, store all the paths
                            pathlist = list()
                            # A full taxonomy list look like this:
                            # [root', 'cellular organisms', 'Bacteria', 'Proteobacteria', 'Gammaproteobacteria',
                            # 'Enterobacterales', 'Enterobacteriaceae', 'Salmonella', 'Salmonella enterica',
                            # 'Salmonella enterica subsp. salamae'] - the first three entries can be ignored
                            for taxlevel in taxonomylist[3:]:
                                # Stop at the species level e.g. ignore Salmonella enterica subsp. salamae
                                if len(str(taxlevel).split(' ')) <= 2:
                                    # Add the taxonomy level to the path e.g. /outputpath/Bacteria/Proteobacteria
                                    level = ''.join(char for char in taxlevel if char.isalnum() or char == ' ')
                                    level = level.replace(' ', '_').replace('.', '')
                                    path = os.path.join(path, level)
                                    # Add this path to the list of paths
                                    pathlist.append(path)
                            # Iterate through all the paths for the accession, create the path if necessary, and create
                            # the .fasta file containing all sequences belonging to this taxonomic level
                            for taxpath in pathlist:
                                make_path(taxpath)
                                #
                                originalpath = os.path.split(taxpath)[-1]
                                name = originalpath + '_{}'.format(self.analysistype)
                                if len(str(originalpath).split('_')) == 2 and 'group' not in originalpath \
                                        and 'unclassified' not in originalpath:
                                    genus = originalpath.split('_')[0]
                                    genusfolder = os.path.join(self.genuspath, genus)
                                    make_path(genusfolder)
                                    genusfile = os.path.join(genusfolder, genus + '.fa')
                                    with open(genusfile, 'a') as genera:
                                        SeqIO.write(self.recorddict[self.linkingdict[accession]], genera, 'fasta')
                                #
                                pathfile = os.path.join(os.path.dirname(taxpath), name)
                                if not os.path.isfile(pathfile):
                                    # Create a metadata object for each file
                                    metadata = MetadataObject()
                                    metadata.originalpath = originalpath
                                    metadata.name = name
                                    metadata.pathfile = pathfile
                                    metadata.childpaths = taxpath.replace(self.output, '').lstrip(os.path.sep)
                                    self.samples.add(metadata)
                                # Open the .fasta file to append
                                with open(pathfile, 'a') as filtered:
                                    # Use SeqIO to write the record stored in recorddict using a key of the accession
                                    # in self.linkingdict
                                    # e.g. accession = NR_044371.1,
                                    # self.linkingdict[accession] = gi|343205886|ref|NR_044371.1|
                                    # recorddict[self.linkingdict[accession]] = record of gi|343205886|ref|NR_044371.1|
                                    SeqIO.write(self.recorddict[self.linkingdict[accession]], filtered, 'fasta')
                                # if len(str(taxlevel).split(' ')) == 2:
                                #     print(metadata.datastore)
                            with open(bacteriafile, 'a') as bacteria:
                                SeqIO.write(self.recorddict[self.linkingdict[accession]], bacteria, 'fasta')
                    except ValueError:
                        pass
                        print('missing', accession, taxid)
            quit()
            # Pickle and write the objects to file
            with open(self.pathobject, 'wb') as pathfile:
                pickle.dump(self.samples, pathfile)
                # Load the pickled objects
        else:
            with open(self.pathobject, 'rb') as pathfile:
                self.samples = pickle.load(pathfile)

    def cluster(self):
        """

        """
        from threading import Thread
        printtime('Clustering sequences to remove redundancies', self.starttime)
        make_path(self.clusterpath)
        if not os.path.isfile(self.clusterobject):
            # Create the threads for the analysis
            for _ in range(self.cpus):
                threads = Thread(target=self.clusterthreads, args=())
                threads.setDaemon(True)
                threads.start()

            for sample in self.samples:
                sample.cluster = GenObject()
                # print(sample.datastore)
                sample.cluster.outputpath = os.path.join(self.clusterpath, os.path.dirname(sample.childpaths))
                make_path(sample.cluster.outputpath)
                sample.cluster.outputfile = os.path.join(sample.cluster.outputpath, sample.name)
                sample.cluster.clusterfile = sample.cluster.outputpath + '.clstr'
                sample.cluster.call = ['cd-hit-est', '-i', sample.pathfile, '-c', '{}'.format(self.cutoff * 0.01), '-o',
                                       sample.cluster.outputfile]
                self.queue.put(sample)
            self.queue.join()
            # Pickle and write the objects to file
            with open(self.clusterobject, 'wb') as cluster:
                pickle.dump(self.samples, cluster)

        else:
            # Load the pickled objects
            with open(self.clusterobject, 'rb') as cluster:
                self.samples = pickle.load(cluster)

    def clusterthreads(self):
        import subprocess
        while True:
            sample = self.queue.get()
            if not os.path.isfile(sample.cluster.outputfile):
                #
                subprocess.call(sample.cluster.call, stdout=self.devnull, stderr=self.devnull)
            self.queue.task_done()

    def clustertargets(self):
        """

        """
        printtime('Creating taxonomically-sorted targets', self.starttime)
        # Remove the folder due to the appending nature of the file creation
        try:
            shutil.rmtree(self.targetpath)
            # shutil.rmtree(self.baitpath)
        # If the folder doesn't exist, ignore
        except FileNotFoundError:
            pass
        make_path(self.targetpath)
        # make_path(self.baitpath)
        # for sample in self.samples:
        #     if len(sample.originalpath.split('_')) == 2:
        #         print(sample.datastore)
        #         with open(sample.pathfile, 'r') as test:
        #             for record in test:
        #                 print(record)
        #     quit()
        # quit()
        for sample in self.samples:
            sample.targetpath = os.path.join(self.targetpath, os.path.dirname(sample.childpaths))
            make_path(sample.targetpath)
            sample.targetfile = os.path.join(sample.targetpath,
                                             '{:03d}_{}_probe.fa'
                                             .format(len(sample.childpaths.split(os.path.sep)), self.analysistype))

            # baitfile = os.path.join(self.baitpath, '{}_bait.fa'.format(self.analysistype))
            # with open(baitfile, 'a') as bait:
            with open(sample.targetfile, 'a') as combined:
                if len(sample.originalpath.split('_')) == 2:
                    record = list(SeqIO.parse(sample.pathfile, 'fasta'))
                    SeqIO.write(record, combined, 'fasta')
                    # SeqIO.write(record, bait, 'fasta')
                else:
                    with open(sample.cluster.outputfile, 'r') as clusters:
                        count = 1
                        for record in SeqIO.parse(clusters, 'fasta'):
                            record.id = '{}_{:02d}'.format(sample.name, count)
                            record.name = ''
                            record.description = ''
                            record.seq = record.seq
                            SeqIO.write(record, combined, 'fasta')
                            # SeqIO.write(record, bait, 'fasta')
                            count += 1

    def align(self):
        """
        Use clustal omega in BioPython to align the filtered, sorted multifasta files
        """
        # from Bio.Align.Applications import ClustalOmegaCommandline
        from Bio.Align.Applications import MuscleCommandline
        from threading import Thread
        printtime('Aligning alleles', self.starttime)
        if not os.path.isfile(self.alignobject):
            # Create the threads for the analysis
            for _ in range(self.cpus):
                threads = Thread(target=self.alignthreads, args=())
                threads.setDaemon(True)
                threads.start()
            for sample in self.samples:
                sample.alignpath = os.path.join(self.path, 'alignedalleles', '')
                make_path(sample.alignpath)
                # Create a list to store objects
                sample.alignedalleles = list()
                sample.aligned = os.path.join(sample.alignpath, sample.filename)
                sample.alignedalleles.append(sample.aligned)
                #
                muscle = MuscleCommandline(input=sample.pathfile,
                                           out=sample.aligned,
                                           maxiters=5)
                # Save a string of the command to the object
                sample.muscle = str(muscle)
                self.queue.put((sample, muscle, sample.pathfile, sample.aligned))
            self.queue.join()
            # Pickle and write the objects to file
            with open(self.alignobject, 'wb') as alignfile:
                pickle.dump(self.samples, alignfile)
        else:
            # Load the pickled objects
            with open(self.alignobject, 'rb') as alignfile:
                self.samples = pickle.load(alignfile)

    def alignthreads(self):
        while True:
            sample, clustalomega, outputfile, aligned = self.queue.get()
            if not os.path.isfile(aligned):
                # Perform the alignments
                # noinspection PyBroadException
                try:
                    clustalomega()
                # Files with a single sequence cannot be aligned. Copy the original file over to the aligned folder
                except Exception:
                    shutil.copyfile(outputfile, aligned)
            self.queue.task_done()

    def probefinder(self):
        """
        Find the longest probe sequences
        """
        from Bio import AlignIO
        from Bio.Align import AlignInfo
        import numpy
        printtime('Finding and filtering probe sequences', self.starttime)
        if not os.path.isfile(self.probeobject):
            for sample in self.samples:
                # A list to store the metadata object for each alignment
                sample.gene = list()
                for align in sample.alignedalleles:
                    # Create an object to store all the information for each alignment file
                    metadata = GenObject()
                    metadata.name = os.path.basename(align).split('.')[0]
                    metadata.alignmentfile = align
                    # Create an alignment object from the alignment file
                    metadata.alignment = AlignIO.read(align, 'fasta')
                    metadata.summaryalign = AlignInfo.SummaryInfo(metadata.alignment)
                    # The dumb consensus is a very simple consensus sequence calculated from the alignment. Default
                    # parameters of threshold=.7, and ambiguous='X' are used
                    consensus = metadata.summaryalign.dumb_consensus()
                    metadata.consensus = str(consensus)
                    # The position-specific scoring matrix (PSSM) stores the frequency of each based observed at each
                    # location along the entire consensus sequence
                    metadata.pssm = metadata.summaryalign.pos_specific_score_matrix(consensus)
                    metadata.identity = list()
                    # Find the prevalence of each base for every location along the sequence
                    for line in metadata.pssm:
                        # Calculate the frequency of the most common base - only count A, C, G, and T
                        try:
                            bases = [line['A'], line['C'], line['G'], line['T']]
                        except KeyError:
                            bases = [line['a'], line['c'], line['g'], line['t']]
                        metadata.identity.append(
                            float('{:.2f}'.format(max(bases) / sum([value for value in line.values()]) * 100)))
                    # List to store metadata objects
                    metadata.windows = list()
                    # Variable to store whether a suitable probe has been found for the current organism + gene pair.
                    # As the probe sizes are evaluated in descending size, as soon as a probe has been discovered, the
                    # search for more probes can stop, and subsequent probes will be smaller than the current one(s)
                    passing = False
                    # Create sliding windows of size self.max - self.min from the list of identities for each column
                    # of the alignment
                    for i in reversed(range(self.min, self.max + 1)):
                        if not passing:
                            windowdata = MetadataObject()
                            windowdata.size = i
                            windowdata.max = 0
                            windowdata.sliding = list()
                            # Create a counter to store the starting location of the window in the sequence
                            n = 0
                            # Create sliding windows from the range of sizes for the list of identities
                            windows = self.window(metadata.identity, i)
                            # Go through each window from the collection of sliding windows to determine which window(s)
                            # has (have) the best results
                            for window in windows:
                                # Create another object to store all the data for the window
                                slidingdata = MetadataObject()
                                # Only consider the window if every position has a % identity greater than the cutoff
                                if min(window) > self.cutoff:
                                    # Populate the object with the necessary variables
                                    slidingdata.location = '{}:{}'.format(n, n + i)
                                    slidingdata.min = min(window)
                                    slidingdata.mean = float('{:.2f}'.format(numpy.mean(window)))
                                    slidingdata.sequence = str(consensus[n:n+i])
                                    # Create attributes for evaluating windows. A greater/lesser
                                    # windowdata.max/windowdata.min means better/lesser overall % identity, respectively
                                    windowdata.max = slidingdata.mean if slidingdata.mean >= windowdata.max \
                                        else windowdata.max
                                    windowdata.min = slidingdata.mean if slidingdata.mean <= windowdata.max \
                                        else windowdata.min
                                    # Add the object to the list of objects
                                    windowdata.sliding.append(slidingdata)
                                    passing = True
                                n += 1
                            # All the object to the list of objects
                            metadata.windows.append(windowdata)
                    # All the object to the list of objects
                    sample.gene.append(metadata)
            # Pickle and write the objects to file
            with open(self.probeobject, 'wb') as probefile:
                pickle.dump(self.samples, probefile)
        else:
            # Load the pickled objects
            with open(self.probeobject, 'rb') as probefile:
                self.samples = pickle.load(probefile)

    def probes(self):
        """
        Find the 'best' probes for each gene by evaluating the percent identity of the probe to the best recorded
        percent identity for that organism + gene pair
        """
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        from Bio.SeqRecord import SeqRecord
        printtime('Determining optimal probe sequences', self.starttime)
        make_path(self.probepath)
        # Remove the combined file
        try:
            os.remove(os.path.join(self.probepath, 'gdcs_combined.fasta'))
        except FileNotFoundError:
            pass
        for sample in self.samples:
            # Make a folder to store the probes
            sample.gdcsoutputpath = self.probepath
            sample.gdcscombined = os.path.join(sample.gdcsoutputpath, 'gdcs_combined.fasta')
            with open(sample.gdcscombined, 'a') as combined:
                for gene in sample.gene:
                    # Open the file to append
                    gene.gdcsoutputfile = os.path.join(sample.gdcsoutputpath, '{}_gdcs.tfa'.format(gene.name))

                    for window in gene.windows:
                        # Variable to record whether a probe has already been identified from this gene
                        passed = False
                        for sliding in window.sliding:
                            # Only consider the sequence if the sliding object has data, if the probe in question
                            # has a mean identity equal to the highest observed identity for that probe size, and
                            # if the mean identity is greater or equal than the lowest observed identity
                            if sliding.datastore and sliding.mean == window.max and sliding.mean >= window.min \
                                    and not passed:
                                dnaseq = Seq(sliding.sequence, IUPAC.unambiguous_dna)
                                # Create a sequence record using BioPython
                                gene.fasta = SeqRecord(dnaseq,
                                                       #  Use the gene name as the header
                                                       id=gene.name,
                                                       description='')
                                # Write each probe to the files
                                if not os.path.isfile(gene.gdcsoutputfile):
                                    with open(gene.gdcsoutputfile, 'w') as allelefile:
                                        SeqIO.write(gene.fasta, allelefile, 'fasta')
                                SeqIO.write(gene.fasta, combined, 'fasta')
                                passed = True

    def targets(self):
        """

        """
        printtime('Creating taxonomically-sorted targets', self.starttime)
        # Remove the folder due to the appending nature of the file creation
        try:
            shutil.rmtree(self.targetpath)
        # If the folder doesn't exist, ignore
        except FileNotFoundError:
            pass
        make_path(self.targetpath)
        #
        for sample in self.samples:
            sample.targetpath = os.path.join(self.targetpath, os.path.dirname(sample.childpaths))
            make_path(sample.targetpath)
            sample.targetfile = os.path.join(sample.targetpath,
                                             '{:03d}_{}_probe.fa'
                                             .format(len(sample.childpaths.split(os.path.sep)), self.analysistype))
            #
            for gene in sample.gene:
                try:
                    with open(sample.targetfile, 'a') as targetfile:
                        SeqIO.write(gene.fasta, targetfile, 'fasta')
                except KeyError:
                    print(sample.name, sample.targetpath, sample.targetfile)

    @staticmethod
    def window(iterable, size):
        """
        https://coderwall.com/p/zvuvmg/sliding-window-in-python
        :param iterable: string from which sliding windows are to be created
        :param size: size of sliding window to create
        """
        i = iter(iterable)
        win = []
        for e in range(0, size):
            win.append(next(i))
        yield win
        for e in i:
            win = win[1:] + [e]
            yield win

    def __init__(self, args):
        """
        :param args: command line arguments
        """
        import multiprocessing
        from queue import Queue
        # Initialise variables
        self.starttime = args.start
        self.min = args.min
        self.max = args.max
        self.cutoff = int(args.cutoff)
        self.analysistype = args.analysis
        # Define variables based on supplied arguments
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), 'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.targetfile = os.path.join(self.path, args.targetfile)
        #
        assert os.path.isfile(self.targetfile), 'Target file does not exist in the path {0!r:s}' \
            .format(self.targetfile)
        self.output = os.path.join(self.path, 'filtered')
        self.taxidfile = os.path.join(self.path, args.database)
        self.reducedtaxidfile = self.taxidfile + '_reduced'
        self.completefile = self.path + 'complete.txt'
        self.objectpath = os.path.join(self.path, 'objects')
        self.pathobject = os.path.join(self.objectpath, 'pathobject')
        make_path(self.objectpath)
        self.alignobject = os.path.join(self.objectpath, 'alignobject')
        self.probeobject = os.path.join(self.objectpath, 'probeobject')
        self.clusterpath = os.path.join(self.path, 'clustered')
        self.clusterobject = os.path.join(self.objectpath, 'clusterobject')
        self.targetpath = os.path.join(self.path, 'targets')
        self.baitpath = os.path.join(self.path, 'bait')
        self.genuspath = os.path.join(self.path, 'genera')
        assert os.path.isfile(self.taxidfile) or os.path.isfile(self.reducedtaxidfile), \
            'Cannot find nucl_gb.accession2taxid file in the path. Please ensure that you ' \
            'downloaded, decompressed, and placed into the supplied path the file from ' \
            'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        self.probepath = os.path.join(self.path, 'probes')
        # Dictionary to store indexed SeqIO records
        self.recorddict = dict()
        # Dictionary to link :self.recorddict to accessions
        self.linkingdict = dict()
        self.samples = set()
        self.cpus = multiprocessing.cpu_count()
        self.queue = Queue(maxsize=self.cpus)
        self.devnull = open(os.devnull, 'w')
        # Filter the input file
        self.runner()

if __name__ == '__main__':
    import time
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='')
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
    parser.add_argument('-m', '--min',
                        default=20,
                        help='Minimum size of probe to create')
    parser.add_argument('-M', '--max',
                        default=100,
                        help='Maximum size of probe to create')
    parser.add_argument('-c', '--cutoff',
                        default=95,
                        help='Cutoff percent identity of a nucleotide location to use')
    parser.add_argument('-a', '--analysis',
                        default='sixteens',
                        help='Specify the gene family: sixteens or cpn60')
    parser.add_argument('-d', '--database',
                        default='nucl_gb.accession2taxid',
                        help='')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.start = time.time()

    # Run the script
    Filter(arguments)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.start) + '\033[0m')

'''
/nas0/bio_requests/8878/cpn60
-t
cpn60_ref_nut_seq.txt
-a
cpn60
-c
95
'''

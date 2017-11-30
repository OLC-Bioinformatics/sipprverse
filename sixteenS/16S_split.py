#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import *
import subprocess
__author__ = 'adamkoziol'


class Split(object):

    def split(self):
        from Bio import SeqIO
        """
        
        """
        printtime('Splitting multifasta file into individual files', self.start)
        with open(self.targetfile, 'r') as target:
            for record in SeqIO.parse(target, 'fasta'):
                sample = MetadataObject()

                sample.name = os.path.splitext(record.id.split('|')[-2])[0]
                sample.outputfile = os.path.join(self.splitpath, sample.name)
                if not os.path.isfile(sample.outputfile):
                    with open(sample.outputfile, 'w') as outputfile:
                        SeqIO.write(record, outputfile, 'fasta')
                self.samples.append(sample)
        self.mashsketch()

    def mashsketch(self):
        """

        """
        from threading import Thread
        printtime('Sketching individual files with mash', self.start)
        # Create the threads for the analysis
        for _ in range(self.cpus):
            threads = Thread(target=self.mashthreads, args=())
            threads.setDaemon(True)
            threads.start()
        filelist = list()
        for sample in self.samples:
            sample.mashoutputfilenoext = os.path.join(self.sketchpath, sample.name)
            sample.mashoutputfile = sample.mashoutputfilenoext + '.msh'
            sample.mashsketch = ['mash', 'sketch',
                                 '-k', '19',
                                 '-s', '100',
                                 sample.outputfile,
                                 '-o', sample.mashoutputfilenoext]
            filelist.append(sample.mashoutputfile)
            # print(sample.name, sample.datastore)
            self.queue.put(sample)
        self.queue.join()
        with open(self.filelist, 'w') as files:
            files.write('\n'.join(filelist))
        self.mashpaste()

    def mashthreads(self):
        while True:
            sample = self.queue.get()
            if not os.path.isfile(sample.mashoutputfile):
                subprocess.call(sample.mashsketch, stdout=self.devnull, stderr=self.devnull)
            self.queue.task_done()

    def mashpaste(self):
        """

        """
        printtime('Merging individual sketch files', self.start)
        # 'cd', self.outputpath, '&&',
        outputfilenoext = os.path.join(self.outputpath, 'merged')
        outputfile = outputfilenoext + '.msh'
        gzipfile = outputfile + '.gz'
        pastecall = ['mash', 'paste', outputfilenoext, '-l', self.filelist]
        if not os.path.isfile(outputfile) and not os.path.isfile(gzipfile):
            #
            subprocess.call(pastecall, stdout=self.devnull, stderr=self.devnull)

    def __init__(self, args):
        """
        :param args: command line arguments
        """
        import multiprocessing
        from queue import Queue
        # Initialise variables
        self.start = args.start
        # Define variables based on supplied arguments
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), 'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.targetfile = os.path.join(self.path, args.targetfile)
        #
        assert os.path.isfile(self.targetfile), 'Target file does not exist in the path {0!r:s}' \
            .format(self.targetfile)
        self.splitpath = os.path.join(self.path, 'split')
        make_path(self.splitpath)
        self.sketchpath = os.path.join(self.path, 'sketches')
        make_path(self.sketchpath)
        self.outputpath = os.path.join(self.path, 'merged')
        make_path(self.outputpath)
        self.samples = list()
        self.filelist = self.path + 'sketches.txt'
        self.cpus = multiprocessing.cpu_count()
        self.queue = Queue(maxsize=self.cpus)
        self.devnull = open(os.devnull, 'w')
        # Filter the input file
        self.split()

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
                        help='Name of sixteenS target file to process. Note that this file must be within the '
                             'supplied path.')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.start = time.time()

    # Run the script
    Split(arguments)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.start) + '\033[0m')

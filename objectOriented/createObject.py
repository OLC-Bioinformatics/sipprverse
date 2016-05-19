#!/usr/bin/env python
from SPAdesPipeline.OLCspades.accessoryFunctions import *

__author__ = 'adamkoziol'


class ObjectCreation(object):
    def createobject(self):
        from glob import glob
        from shutil import copyfile
        # Grab any .fastq files in the path
        fastqfiles = glob('{}*.fastq*'.format(self.path))
        # Extract the base name of the globbed name + path provided
        fastqnames = map(lambda x: os.path.split(x)[1], filer(fastqfiles))
        # Iterate through the names of the fastq files
        for fastqname in sorted(fastqnames):
            # Set the name
            metadata = MetadataObject()
            metadata.name = fastqname
            # Set the destination folder
            outputdir = '{}{}'.format(self.path, fastqname)
            # Make the destination folder
            make_path(outputdir)
            # Get the fastq files specific to the fastqname
            specificfastq = glob('{}{}*.fastq*'.format(self.path, fastqname))
            # Copy the files to the output folder
            try:
                # Copy the .gz files to :self.path/:filename
                for fastqfile in specificfastq:
                    destinationfastq = '{}/{}'.format(outputdir, os.path.split(fastqfile)[1])
                    if not os.path.isfile(destinationfastq):
                        copyfile(fastqfile, destinationfastq)
            # Except os errors
            except OSError as exception:
                # If there is an exception other than the file exists, raise it
                if exception.errno != errno.EEXIST:
                    raise
            # Initialise the general and run categories
            metadata.general = GenObject()
            metadata.run = GenObject()
            # Populate the .fastqfiles category of :self.metadata
            metadata.general.fastqfiles = [fastq for fastq in glob('{}/{}*.fastq*'.format(outputdir, fastqname))
                                           if 'trimmed' not in fastq]
            # Add the output directory to the metadata
            metadata.general.outputdirectory = outputdir
            # Append the metadata to the list of samples
            self.samples.append(metadata)

    def __init__(self, inputobject):
        self.samples = []
        self.path = inputobject.sequencepath
        self.createobject()

#!/usr/bin/env python
from SPAdesPipeline.OLCspades.fastqCreator import *

__author__ = 'adamkoziol'


class FastqCreate(CreateFastq):
    def fastqmover(self):
        """Links .fastq files created above to :sequencepath"""
        from re import sub
        import errno
        # Create the project path variable
        self.projectpath = self.fastqdestination + "/Project_" + self.projectname
        # Create the sequence path if necessary
        make_path(self.sequencepath)
        # Iterate through all the sample names
        for sample in self.metadata.samples:
            # Make the outputdir variable
            outputdir = '{}{}'.format(self.sequencepath, sample.name)
            # Glob all the .gz files in the subfolders - projectpath/Sample_:sample.name/*.gz
            for fastq in sorted(glob('{}/Sample_{}/*.gz'.format(self.projectpath, sample.name))):
                # Try/except loop link .gz files to self.path
                try:
                    # Symlink fastq file to the path, but renames them first using the sample number.
                    # 2015-SEQ-1283_GGACTCCT-GCGTAAGA_L001_R1_001.fastq.gz is renamed:
                    # 2015-SEQ-1283_S1_L001_R1_001.fastq.gz
                    make_path(outputdir)
                    os.symlink(fastq, '{}/{}'.format(outputdir, os.path.basename(sub('\w{8}-\w{8}',
                                                                                     'S{}'.format(
                                                                                         sample.run.SampleNumber),
                                                                                     fastq))))
                # Except os errors
                except OSError as exception:
                    # If there is an exception other than the file exists, raise it
                    if exception.errno != errno.EEXIST:
                        raise
            # Repopulate .strainfastqfiles with the freshly-linked files
            fastqfiles = glob('{}/{}*.fastq*'.format(outputdir, sample.name))
            fastqfiles = [fastq for fastq in fastqfiles if 'trimmed' not in fastq]
            # Populate the metadata object with the name/path of the fastq files
            sample.general.fastqfiles = fastqfiles
            # Save the outputdir to the metadata object
            sample.run.outputdirectory = outputdir
            sample.general.bestassemblyfile = True
            sample.general.trimmedcorrectedfastqfiles = sample.general.fastqfiles
            sample.commands = GenObject()

    def __init__(self, inputobject):
        self.sequencepath = inputobject.sequencepath
        CreateFastq.__init__(self, inputobject)

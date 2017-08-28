#!/usr/bin/env python
import sys

import Bio

from SPAdesPipeline.OLCspades.accessoryFunctions import *
from SPAdesPipeline.OLCspades.bowtie import *

__author__ = 'adamkoziol'


class Versions(object):
    def versions(self):
        for sample in self.metadata:
            # Initialise the attribute
            sample.software = GenObject()
            # Populate the versions of the software used
            ss = sample.software
            ss.python = self.python
            ss.arch = self.arch
            ss.bowtie2 = self.bowversion
            ss.samtools = self.samversion
            ss.biopython = self.biopython
            ss.mira = self.mira

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.commit = inputobject.commit
        # Determine the versions of the software used
        printtime('Populating metadata', self.start)
        self.python = sys.version.replace('\n', '')
        self.arch = ", ".join(os.uname())
        self.bowversion = Bowtie2CommandLine(version=True)()[0].split('\n')[0].split()[-1]
        self.samversion = get_version(['samtools', '--version']).split('\n')[0].split()[1]
        self.biopython = Bio.__version__
        self.mira = get_version(['mirabait']).split('\n')[0].split()[3].replace(')', '')
        self.versions()

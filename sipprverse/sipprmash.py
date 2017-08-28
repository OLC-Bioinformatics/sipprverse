#!/usr/bin/env python
from SPAdesPipeline.OLCspades.mash import *

__author__ = 'adamkoziol'


class SipprMash(Mash):
    def __init__(self, inputobject, analysistype):
        self.runmetadata = inputobject.runmetadata
        # ref file path is used to work with sub module code with a different naming scheme
        self.reffilepath = inputobject.targetpath
        self.reportpath = inputobject.reportpath
        self.starttime = inputobject.starttime
        self.cpus = inputobject.cpus
        Mash.__init__(self, self, analysistype)

#!/usr/bin/env python
from sipprCommon.sippingmethods import *
__author__ = 'adamkoziol'


class Custom(object):

    def reporter(self):
        """
        Calls the sipping method and creates a report of the results
        """
        #
        Sippr(self)
        #
        header = 'Strain,Gene,PercentIdentity,FoldCoverage\n'
        data = ''
        with open('{}/{}.csv'.format(self.reportpath, self.analysistype), 'w') as report:
            for sample in self.metadata:
                data += sample.name + ','
                multiple = False
                for name, identity in sample[self.analysistype].results.items():
                    if not multiple:
                        data += '{},{},{}\n'.format(name, identity.items()[0][0], identity.items()[0][1])
                    else:
                        data += ',{},{},{}\n'.format(name, identity.items()[0][0], identity.items()[0][1])
                    multiple = True
            report.write(header)
            report.write(data)

    # noinspection PyDefaultArgument
    def __init__(self, inputobject, analysistype, cutoff=0.98, matchbonus=2, builddict=dict(), extension='.bt2'):
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        try:
            self.targetpath = inputobject.customtargetpath if inputobject.customtargetpath else inputobject.targetpath
        except AttributeError:
            self.targetpath = inputobject.targetpath
        self.reportpath = os.path.join(self.path, 'reports')
        make_path(self.reportpath)
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.analysistype = analysistype
        self.cpus = inputobject.cpus
        self.homepath = inputobject.homepath
        self.cutoff = cutoff
        self.matchbonus = matchbonus
        self.builddict = builddict
        self.bowtiebuildextension = extension

#!/usr/bin/env python
from SPAdesPipeline.OLCspades.accessoryFunctions import *

__author__ = 'adamkoziol'


class VtypeResults(object):
    def reporter(self):
        printtime('Creating reports', self.start)
        # Create the report path if necessary
        make_path(self.reportpath)
        # Initialise strings to hold results as well as file names
        detailedheader = 'Strain,MatchName,Vtype,PercentIdentity,AverageCoverage,SNPs,Gaps\n'
        detailedresults = ''
        detailedresultsfile = '{}/{}_detailedResults.csv'.format(self.reportpath, self.analysistype)
        header = 'Strain,virulenceType\n'
        results = ''
        resultsfile = '{}/{}_results.csv'.format(self.reportpath, self.analysistype)
        # Populate the strings for each sample
        for sample in self.metadata:
            # Initialise the set containing the virulence types
            virulenceset = set()
            # Ensure that there are results in the sample object
            if 'results' in sample[self.analysistype].datastore:
                if sample[self.analysistype].results != 'NA':
                    # Iterate through the results dictionary
                    for vtx in sample[self.analysistype].results:
                        # Create easier to understand variable names
                        percentidentity = sample[self.analysistype].results[vtx].items()[0][0]
                        averagedepth = sample[self.analysistype].results[vtx].items()[0][1]
                        snps = str(sample[self.analysistype].resultssnp[vtx])
                        gaps = str(sample[self.analysistype].resultsgap[vtx])
                        # print sample.name, vtx, percentidentity, averagedepth, snps, gaps
                        # The alleles have names like stx2a_2_F_4, this splicing yields stx2a
                        allelename = vtx[:5] + "_F" if "F" in vtx else vtx[:5] + "_R"
                        # Add the alleles to a set
                        virulenceset.add(allelename)
                        # Populate the detailed results with the appropriate values
                        detailedresults += ','.join([sample.name, vtx, allelename, percentidentity,
                                                     averagedepth, snps, gaps])
                        # Start a new line for the next set of results
                        detailedresults += '\n'
                    # Add an extra line between samples
                    detailedresults += '\n'
                # Initialise the set to hold the vtypes
                vtxset = set()
                for vtx in virulenceset:
                    vtype = vtx.split('_')[0]
                    # If forward and reverse primers are present (eg stx2c_F and stx2c_R), then add the vtype to the set
                    if vtype + '_F' in virulenceset and vtype + '_R' in virulenceset:
                        vtxset.add(vtype)
                # Populate the results
                if vtxset:
                    results += sample.name + ','
                    results += ';'.join(sorted(vtxset))
                    results += '\n'
                else:
                    results += sample.name + ','
                    results += '-\n'
                # Set attributes for metadata collection
                sample[self.analysistype].detailedresultsfile = detailedresultsfile
                sample[self.analysistype].resultsfile = resultsfile
                sample[self.analysistype].reportdir = self.reportpath
                sample.general.vtxset = sorted(vtxset)
                sample.general.vtype = ','.join(sorted(vtxset))
        # Open and write the detail and regular reports
        with open(detailedresultsfile, 'wb') as writedetails:
            writedetails.write(detailedheader)
            writedetails.write(detailedresults)
        with open(resultsfile, 'wb') as writeresults:
            writeresults.write(header)
            writeresults.write(results)

    def __init__(self, inputobject, analysistype):
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        self.targetpath = inputobject.customtargetpath
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.analysistype = analysistype
        self.reportpath = os.path.join(self.path, 'reports')
        # Create the reports
        self.reporter()

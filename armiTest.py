from glob import glob
from collections import defaultdict
from Bio import SeqIO
import os, json, pysam, pysamstats, re, time, sys


# Initialise the count used in the dotter function
count = 0


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    # Use a global variable
    global count
    if count <= 80:
        sys.stdout.write('.')
        count += 1
    else:
        sys.stdout.write('\n[%s] .' % (time.strftime("%H:%M:%S")))
        count = 1


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

holdingDict = defaultdict(make_dict)
seqDict = defaultdict(make_dict)
plusdict = defaultdict(make_dict)
start = time.time()

# targets = glob("/media/nas/akoziol/Pipeline_development/GeneSipperV2/baitTest/targets/ARMI/*.fa")
strain = "1042-New"
bamFile = "/media/nas/akoziol/Pipeline_development/GeneSipperV2/ARMI/ARMI_sorted.bam"
target = "/media/nas/akoziol/Pipeline_development/GeneSipperV2/ARMI/ARMIBait.fa"
splitTargets = glob("/media/nas/akoziol/Pipeline_development/GeneSipperV2/baitTest/targets/ARMI/*.fa")
JSONProfile = "/media/nas/akoziol/Pipeline_development/GeneSipperV2/ARMI/ARMIparse.json"
geneProfile = "/media/nas/akoziol/Pipeline_development/GeneSipperV2/ARMI/targetProfile.json"


def populateDict():
    global holdingDict
    if not os.path.isfile(JSONProfile):
        for rec in pysamstats.stat_variation(alignmentfile=bamFile, fafile=target):
            # holdingDict[rec['chrom'].split("_")[0]] = {}

            # geneName = rec['chrom'].split(".")[0]
            aros = re.findall("(?<=ARO:)3\d{6}", rec['chrom'])
            aro = " ".join(aros)
            # if len(aros) > 1:
            #     print str(aros)
            # aros = rec['chrom'].split("!")[1].split("_")[0]
            if rec['matches'] - rec['mismatches'] > 0:
                holdingDict[rec['chrom']][int(rec['pos'])] = rec['reads_all']
            dotter()


        JSONreport = open(JSONProfile, "wb")
        output = json.dumps(holdingDict, sort_keys=True, indent=4, separators=(',', ': '))
        JSONreport.write(output)
        JSONreport.close()
    else:
        with open(JSONProfile, "rb") as jsonReport:
            # Load the data
            # print strain, targetName, JSONProfile, jsonReport
            holdingDict.update(json.load(jsonReport))


def strainTargetSequenceThreads():
    # print json.dumps(seqDict, sort_keys=True, indent=4, separators=(',', ': '))
    global splitTargets, geneProfile, seqDict
    if not os.path.isfile(geneProfile):
        for target in splitTargets:
            targetName = os.path.basename(target).split(".")[0]
            # print strain, target
            # print strain, target
            handle = open(target, "rU")
            for record in SeqIO.parse(handle, "fasta"):
                # print strain, target, record.id, len(record)
                # print type(record.id), type(record.seq), record.seq
                seqDict[record.id.split(" ")[0]][targetName] = len(record.seq)
                # profileDict[record.id] = str(record.seq)
        JSONreport = open(geneProfile, "wb")
        output = json.dumps(seqDict, sort_keys=True, indent=4, separators=(',', ': '))
        JSONreport.write(output)
        JSONreport.close()
    else:
        with open(geneProfile, "rb") as jsonReport:
            # Load the data
            # print strain, targetName, JSONProfile, jsonReport
            seqDict.update(json.load(jsonReport))


import operator
from ARMICARD import decipher

# targetName = os.path.basename(target).split(".")[0]

def parseDict():
    global holdingDict
    global seqDict
    global plusdict
    global strain
    for gene in seqDict:
        for aros in seqDict[gene]:
            length = seqDict[gene][aros]
            # print gene, length, aros
            # aro = aros.split(" ")
            # for ar in aro:
            #     # pass
            #     # if len(aro) > 1:
            matches = 0
            totalDepth = 0
            minDepth = 10
            # for presence in sorted(holdingDict[gene].items(), key=operator.itemgetter(0)):
            for presence in sorted(holdingDict[gene]):
                depth = holdingDict[gene][presence]
                matches += 1
                totalDepth += depth
                if depth < minDepth:
                    minDepth = depth
                # if gene == "AP009048.1.gene3309" and aros == "3000502":
                #     print gene, aros, length, presence, type(presence), matches, totalDepth
            averageDepth = float("%0.2f" % (float(totalDepth) / float(length)))
            percentID = float("%0.2f" % (float(matches) / float(length))) * 100
            # print gene, aros, matches, length, percentID, averageDepth
            if percentID > 70 and minDepth > 4:
                plusdict[strain][aros] = ["+"]
                # print gene, aros, matches, length, percentID, averageDepth
            else:
                plusdict[strain][aros] = []
    antidict = json.load(open("/media/nas0/Jackson/ARMI_Docker/ARMI/aro3.json"))
    resDict = decipher(plusdict, antidict, "/media/nas/akoziol/Pipeline_development/GeneSipperV2/baitTest/results/armi70_5")
                # print gene, aros, matches, length, percentID, averageDepth
            #         plusdict[ar] = percentID



populateDict()
strainTargetSequenceThreads()
# print json.dumps(holdingDict, sort_keys=True, indent=4, separators=(',', ': '))
print json.dumps(seqDict, sort_keys=True, indent=4, separators=(',', ': '))
parseDict()
# print json.dumps(plusdict, sort_keys=True, indent=4, separators=(',', ': '))

print "\nElapsed Time: %0.2f seconds" % (time.time() - start)

from Bio import SeqIO
from collections import defaultdict
from glob import glob
from Queue import Queue
from threading import Thread
import threading, json, sys, os, time


# Declare queues, list, and dictionaries
dqueue = Queue()
threadlock = threading.Lock()


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


def strainTargetSequenceThreads(seqDict, analysisType):
    # print json.dumps(seqDict, sort_keys=True, indent=4, separators=(',', ': '))
    for strain in seqDict:
        try:
            allelePath = os.path.split(seqDict[strain]["targets"][analysisType][0])[0]
            JSONalleles = "%s/%s_alleles.json" % (allelePath, analysisType)
            if not os.path.isfile(JSONalleles):
                for target in seqDict[strain]["targets"][analysisType]:
                    targetName = os.path.basename(target).split(".")[0]
                    # print strain, target, targetName, analysisType
                    # print strain, target
                    handle = open(target, "rU")
                    for record in SeqIO.parse(handle, "fasta"):
                        # print strain, target, record.id, len(record)
                        # print type(record.id), type(record.seq), record.seq
                        seqDict[strain]["targetSequences"][analysisType][target]["allele"][record.id] = str(record.seq)
        except AttributeError:
            pass
                    # profileDict[record.id] = str(record.seq)
            # JSONreport = open(JSONalleles, "wb")
            # output = json.dumps(profileDict, sort_keys=True, indent=4, separators=(',', ': '))
            # JSONreport.write(output)
            # JSONreport.close()
        # else:
        #     with open(JSONalleles, "rb") as jsonReport:
        #         Load the data
                # profileData = json.load(jsonReport)
                # print json.dumps(profileData, sort_keys=True, indent=4, separators=(',', ': '))
                # for target in seqDict[strain]["targets"][analysisType]:
                #     targetName = os.path.basename(target).split(".")[0]
                #     seqDict[strain]["targetSequences"][analysisType][target]["allele"] = profileData.keys()
                    # for record in profileData:
                    #     print seqDict[strain]
                        # if targetName == record.split("_")[0]:
                        #     seqDict[strain]["targetSequences"][analysisType][target]["allele"][record] = profileData[record]
                            # print strain, target, record, profileData[record]
    # for strain in seqDict:
    # print "!!!", json.dumps(seqDict, sort_keys=True, indent=4, separators=(',', ': '))
        dotter()
    return seqDict
                # Send the threads to makeblastdb
    #             threads = Thread(target=bamParse, args=(dqueue,))
    #             # Set the daemon to true - something to do with thread management
    #             threads.setDaemon(True)
    #             # Start the threading
    #             threads.start()
    #         handle.close()
    # for strain in seqDict:
    #     for target in seqDict[strain]["targets"][analysisType]:
    #         handle = open(target, "rU")
    #         for record in SeqIO.parse(handle, "fasta"):
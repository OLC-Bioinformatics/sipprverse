__author__ = 'akoziol'

import subprocess, os, sys, time, re, errno, logging
# Multiprocessing module
from multiprocessing import Pool


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


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def bowtie2indexTargetsProcesses(targets, targetPath):
    """Indexing multiprocessing helper function"""
    # print '\nIndexing target files'
    # Initialise the args list
    indexArgs = []
    # Initialise the pool of processes - it defaults to the number of processors
    indexPool = Pool()
    for target in targets:
        indexArgs.append((target, targetPath))
    indexPool.map(bowtie2indexTargets, indexArgs)


def bowtie2indexTargets((target, targetPath)):
    # Format the target names properly
    filename = os.path.split(target)[1]
    fileNoExt = filename.split(".")[0]
    # Index the appropriate files
    if not os.path.isfile("%s/%s.1.bt2" % (targetPath, fileNoExt)):
        # Define the indexing command
        indexCommand = "cd %s && bowtie2-build %s %s" % (targetPath, target, fileNoExt)
        # Define /dev/null
        FNULL = open(os.devnull, 'wb')
        # Run the command
        subprocess.call(indexCommand, shell=True, stdout=FNULL, stderr=FNULL)
    dotter()


def bowtie2mappingProcesses(seqDict, analysisType, mapperName):
    """Mapping threads!"""
    # print '\nPerforming bowtie2 reference mapping'
    mappingProcessesArgs = []
    mappingProcessesPool = Pool()
    # uses kmer, targets, readLength, foldCoverage
    for strain in seqDict:
        baitType = seqDict[strain]["bait"]["fastqFiles"].keys()[0]
        for target in seqDict[strain]["targets"][analysisType]:
            baitedFastqList = seqDict[strain]["bait"]["fastqFiles"][baitType]
            mappingProcessesArgs.append((strain, target, mapperName, baitedFastqList, baitType))

    mappingProcessesPool.map(mapping, mappingProcessesArgs)


def mapping((strain, target, mapperName, baitedFastqList, baitType)):
    """Performs the mapping of the  reads to the targets"""
    forwardFastq, reverseFastq = baitedFastqList
    targetNoExt = target.split(".")[0]
    targetName = os.path.basename(target).split(".")[0]
    fastqPath = os.path.split(forwardFastq)[0]
    outPath = "%s/%s" % (fastqPath, targetName)
    make_path(outPath)
    bamName = "%s_%s_%s" % (strain, targetName, mapperName)
    # print strain, targetName, mapperName, baitType, forwardFastq, reverseFastq
#     filename = os.path.split(target)[1]
#     fileNoExt = filename.split(".")[0]
#     megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, fileNoExt, size)
#     filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
#     newPath = "%s/%s" % (filePath, megaName)
#     make_path(newPath)
#     targetPath = "%s/targets/%s/%s_%s" % (path, fileNoExt, fileNoExt, size)
#     # Map the files to the reference

    # bowtie2 -p 8 --very-fast-local -t -a --score-min 'C,0,-1'  -x BACT000002 -1 S-mishmar-haemek-NALR-FFFM_S9_L001_R1_001.fastq.gz -2 S-mishmar-haemek-NALR-FFFM_S9_L001_R2_001.fastq.gz > rmlst_2
    if not os.path.isfile("%s/%s.sam" % (outPath, bamName)):
        # smaltMap = "smalt map -o %s/%s.bam -f bam -n 24 -x %s %s %s" \
        #            % (outPath, bamName, targetNoExt, forwardFastq, reverseFastq)
    # -k 20 --score-min 'C,0,-1' -a --local
        mapCommand = "bowtie2 -p 24 --no-unal --very-sensitive-local -a -x %s -1 %s -2 %s > %s/%s.sam" % (targetNoExt, forwardFastq, reverseFastq, outPath, bamName)
        # print mapCommand
        FNULL = open(os.devnull, 'wb')
        subprocess.call(mapCommand, shell=True, stdout=FNULL, stderr=FNULL)
            # dotter()
    if not os.path.isfile("%s/%s_reFlagged.sam" % (outPath, bamName)):
        modify_bowtie_sam("%s/%s.sam" % (outPath, bamName), 5)
    if not os.path.isfile("%s/%s.bam" % (outPath, bamName)):
        bamCommand = "samtools view -bS %s/%s_reFlagged.sam > %s/%s.bam" % (outPath, bamName, outPath, bamName)
        FNULL = open(os.devnull, 'wb')
        subprocess.call(bamCommand, shell=True, stdout=FNULL, stderr=FNULL)
    # else:
    dotter()

    
def modify_bowtie_sam(bowtie_sam,max_mismatch):
    # fix sam flags for comprehensive pileup
    raw_bowtie_sam = bowtie_sam.split(".")[0]
    with open(bowtie_sam) as sam, open(raw_bowtie_sam + '_reFlagged.sam', 'w') as sam_mod:
        for line in sam:
            if not line.startswith('@'):
                fields = line.split('\t')
                flag = int(fields[1])
                flag = (flag - 256) if (flag & 256) else flag
                m = re.search("NM:i:(\d+)\s",line)
                if m != None:
                    num_mismatch = m.group(1)
                    if int(num_mismatch) <= int(max_mismatch):
                        sam_mod.write('\t'.join([fields[0], str(flag)] + fields[2:]))
                    else:
                        logging.info('Excluding read from SAM file due to missing NM (num mismatches) field: ' + fields[0])
                        num_mismatch = 0
            else:
                sam_mod.write(line)
    return(bowtie_sam,raw_bowtie_sam + '_reFlagged.sam')
from glob import glob
import os, errno, time, subprocess, shutil, re
from multiprocessing import Pool


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def createFastq(miSeqPath, miSeqFolder, path, project, forwardreads, reversereads, customSampleSheet):
    """finds the most recent run folder on the MiSeq and determines which cycle the run is currently on. Runs bcl2fastq
    to create .fastq files with the number of reads specified in SampleSheet.csv"""
    if not miSeqFolder:
        # Get the folders into a list
        miSeqFolders = glob("%s/*/" % miSeqPath)
        # The folder on interest is the most recent one, so sort the folders and choose the last on in the list
        miSeqFolder = sorted(miSeqFolders)[-1]
        folderName = miSeqFolder.split("/")[-2]
    else:
        # Add the path to the folder name
        folderName = miSeqFolder
        miSeqFolder = miSeqPath + miSeqFolder

    # Get the flow cell id from the folder name - it should be the last part of the name, and will be separated by a "-"
    flowcellid = miSeqFolder.split("_")[1]
    # Initialise the folder in which to place the .fastq sequences once they are created
    seqfolder = path + "sequences"
    folderpath = seqfolder + "/" + folderName
    # Make the folderPath path if necessary
    make_path(folderpath)
    # Initialise the reads list to be used for storing the length of the forward and reverse reads
    reads = []
    # Initialise the variables to store data from the original sample sheet
    investigator = ""
    forwardLength = ""
    reverseLength = ""
    projectName = ""
    # bcl2fastq requires an older version of the sample sheet, this recreates the required version
    # Create the new sample sheet
    with open("%s/%s/SampleSheet_modified.csv" % (seqfolder, folderName), "wb") as modifiedSampleSheet:
        # Write the required headings to the file
        modifiedSampleSheet.write("FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n")
        # Open the current sample sheet
        if customSampleSheet:
            # with open("%s/SampleSheet.csv" % customSampleSheet, "rb") as originalSampleSheet:
            originalSampleSheet = open("%sSampleSheet.csv" % customSampleSheet, "rb")
        else:
            # with open("%s/SampleSheet.csv" % miSeqFolder, "rb") as originalSampleSheet:
            originalSampleSheet = open("%s/SampleSheet.csv" % miSeqFolder, "rb")
        for line in originalSampleSheet:
            # As this is a csv file, data are separated by commas - split on commas
            data = line.split(",")
            if "Investigator Name" in line:
                investigator = data[1].rstrip()
            # Here's a Perl-like solution for reading lines after a regex match
            # Perform the regex
            elif "[Reads]" in line:
                # Now read all sublines going forward
                for subline in originalSampleSheet:
                    # Stop reading once "Settings" is matched
                    if "[Settings]" in subline:
                        break
                    reads.append(subline.rstrip())
                # Grab the number of reads in the forward and reverse direction
                forwardLength = int(reads[0].replace(",", "").rstrip())
                reverseLength = int(reads[1].replace(",", "").rstrip())
            # Get the data for each sample from the sheet
            elif "Sample_ID" in line:
                # Iterate through the lines following "Sample_ID"
                for subline in originalSampleSheet:
                    # Split on the ,
                    subdata = subline.split(",")
                    # Populate
                    sampleID = subdata[0]
                    index = "%s-%s" % (subdata[5], subdata[7])
                    # Not fully tested - attempts to catch sample sheets without supplied Sample_Project and/or
                    # Descriptions
                    # Only run this if the project variable has not been provided
                    if not project:
                        # If the there is a project name in the sample sheet
                        if subdata[8]:
                            projectName = subdata[8].rstrip()
                        # Otherwise supply a NA
                        else:
                            projectName = "NA"
                    else:
                        projectName = project
                    # Same as above look for the data, and provide a NA if it's not provided
                    if subdata[9]:
                        description = subdata[9].rstrip()
                    else:
                        description = "NA"
                    # Write the data to the modified sample sheet
                    modifiedSampleSheet.write("%s,1,%s,no_ref,%s,%s,N,NA,%s,%s\n" % (flowcellid, sampleID, index,
                                                    description, investigator, projectName))
    # Calculate the number of cycles completed
    cycles = glob("%s/Data/Intensities/BaseCalls/L001/C*" % miSeqFolder)
    # As the length of reads can be specified as "full", use the forwardLength and reverseLength variables to populate
    # the forwardreads and reversereads as required
    if forwardreads == "full":
        forwardreads = forwardLength
    else:
        forwardreads = int(forwardreads)
    if reversereads == "full":
        reversereads = reversereads
    else:
        reversereads = int(reversereads)
    # Initialise basemask and nohup variables
    basemask = ""
    nohup = ""
    # As the number of cycles required is the number of forward reads + the index(8) + the second index(8) + 1 (just to be safe)
    # Also set the basemask variable as required
    if forwardreads:
        if reversereads:
            readsNeeded = forwardreads + reversereads
            basemask = "Y%sn*,I8,I8,Y%sn*" % (forwardreads, reversereads)
            nohup = "nohup make -j 16"
        else:
            readsNeeded = forwardreads + 16
            basemask = "Y%sn*,I8,I8,n*" % forwardreads
            nohup = "nohup make -j 16 r1"
    else:
        if reverseLength:
            readsNeeded = forwardLength + reverseLength
            basemask = "Y%sn*,I8,I8,Y%sn*" % (forwardLength, reverseLength)
            nohup = "nohup make -j 16"
        else:
            readsNeeded = forwardLength + 16
            basemask = "Y%sn*,I8,I8,n*" % forwardLength
            nohup = "nohup make -j 16 r1"

    # A while loop to calculate when there are enough cycles to perform the analysis as specified in the arguments
    while len(cycles) < readsNeeded:
        print "Currently at %s cycles. Need to wait until the MiSeq reaches cycles %s" % (len(cycles), readsNeeded)
        time.sleep(30)
        cycles = glob("%s/Data/Intensities/BaseCalls/L001/C*" % miSeqFolder)
    if not os.path.isdir("%s/Project_%s" % (folderpath, projectName)):
        # Call configureBclToFastq.pl
        print "Running bcl2fastq"
        bclcall = "configureBclToFastq.pl --input-dir %s/Data/Intensities/BaseCalls/ " \
                  "--output-dir %s --force --sample-sheet %s/SampleSheet_modified.csv " \
                  "--mismatches 1 --no-eamss --fastq-cluster-count 0 --compression none --use-bases-mask %s" \
                  % (miSeqFolder, folderpath, folderpath, basemask)
        print bclcall
        # Define /dev/null
        FNULL = open(os.devnull, 'wb')
        # Run the commands
        subprocess.call(bclcall, shell=True, stdout=FNULL, stderr=FNULL)
        nohupcall = "cd %s/%s && %s" % (seqfolder, folderName, nohup)
        subprocess.call(nohupcall, shell=True, stdout=FNULL, stderr=FNULL)

    # Find all the files created using a nested list comprehension
    # Firstly, list all the directories using os.listdir (the second part of the comprehension)
    # Then look for .gz files in each of the newly discovered folders
    gzFiles = [[sorted(glob("%s/Project_%s/%s/*.gz" % (folderpath, projectName, f)))]
               for f in sorted(os.listdir("%s/Project_%s" % (folderpath, projectName)))]

    # Initialise lists to hold arguments for multiprocessing modules
    bbdukargs = []
    quakeargs = []
    # Iterate through the nested lists
    for i in range(len(gzFiles)):
        # Second level
        for j in range(len(gzFiles[i])):
            # Third level
            for k in range((len(gzFiles[i][j]))):
                # Filenames are strings of the kth element of the jth element of the ith element of gzfiles
                filename = str(gzFiles[i][j][k])
                # Substitute in the sequence number for the indexes used. Should look something like:
                # ..../Sample_2015-SEQ-1106/2015-SEQ-1106_GTAGAGGA-CTAAGCCT_L001_R1_001.fastq.gz becomes
                # Sample_2015-SEQ-1106/2015-SEQ-1106_S1_L001_R1_001.fastq.gz
                subname = re.sub("\w{8}-\w{8}", "S%s" % str(i + 1), filename)
                # Rename the files appropriately
                os.rename(filename, subname)
            # Treat paired vs. unpaired reads differently
            reverse = ""
            cleanreverse = ""
            if len(gzFiles[i][j]) == 2:
                forward, reverse = gzFiles[i][j]
                # Uses lambda to substitute out the "_S1_L001_" and replace the "_001.fastq.gz" portions of the filename
                cleanforward, cleanreverse = map(lambda x: re.sub("_S\d+_L001", "",
                                                       x.replace("_001.fastq.gz", ".clean.fq")), gzFiles[i][j])
            elif len(gzFiles[i][j]) == 1:
                forward = gzFiles[i][j][0]
                cleanforward = re.sub("_S\d+_L001", "", forward.replace("_001.fastq.gz", ".clean.fq"))
            # Append the necessary arguments to the lists for use in multiprocessing later
            bbdukargs.append((forward, reverse, cleanforward, cleanreverse))
            # Create the name of the counts file to be created by jellyfish
            countsfile = "%s/counts" % os.path.split(forward)[0]
            quakeargs.append((cleanforward, cleanreverse, countsfile, seqfolder))
            # Create fastqFiles.txt to store the names of the fastq files to be corrected
            pathname = os.path.split(cleanforward)[0]
            with open("%s/fastqFiles.txt" % pathname, "wb") as fastqFiles:
                fastqFiles.write("%s\t%s" % (cleanforward, cleanreverse))
    # Run multiprocessed bbduker
    bbdukpool = Pool()
    print "Trimming sequences"
    bbdukpool.map(bbduker, bbdukargs)
    # Run the multiprocessed quake
    quakepool = Pool()
    print "Error correcting sequences"
    quakepool.map(quake, quakeargs)


def bbduker((forward, reverse, cleanforward, cleanreverse)):
    """Runs BBduk in a multiprocessed fashion"""
    # Run the command if the output files don't exist
    if not os.path.isfile(cleanforward):
        # Use BBduk to perform quality and adapter trimming. Example command included below:
        # bbduk.sh -Xmx1g in1=2015-SEQ-1132_S1_L001_R1_001.fastq.gz in2=2015-SEQ-1132_S1_L001_R2_001.fastq.gz
        # out1=2015-SEQ-1132_R1.clean.fq out2=2015-SEQ-1132_R2.clean.fq
        # minlen=100 qtrim=rl trimq=20 ktrim=r k=25 mink=11
        # ref=/home/blais/Bioinformatics/bbmap/resources/adapters.fa hdist=1
        if reverse:
            # minlen=100
            bbdukcall = "bbduk.sh -Xmx1g in1=%s in2=%s out1=%s out2=%s qtrim=rl trimq=20 ktrim=r " \
                    "k=25 mink=11 ref=/home/blais/Bioinformatics/bbmap/resources/adapters.fa hdist=1 tpe tbo" \
                    % (forward, reverse, cleanforward, cleanreverse)
        else:
            # minlen=100
            bbdukcall = "bbduk.sh -Xmx1g in=%s out=%s qtrim=rl trimq=20 ktrim=r k=25 mink=11 " \
                        "ref=/home/blais/Bioinformatics/bbmap/resources/adapters.fa hdist=1" % (forward, cleanforward)
        # Define /dev/null
        FNULL = open(os.devnull, 'wb')
        # Run the command - using os.system right now, as subprocess.call failed the last time I tried it
        os.system(bbdukcall)
        # subprocess.call(bbdukcall, shell=True, stdout=FNULL, stderr=FNULL)


def quake((cleanforward, cleanreverse, countsfile, seqFolder)):
    """Runs jellyfish count with the supplied variables in a multiprocessed fashion"""
    correverse = ""
    if cleanreverse:
        # Define the names of the corrected fastq files
        corforward, correverse = map(lambda x: x.replace(".clean.fq", ".clean.cor.fq"), [cleanforward, cleanreverse])
    else:
        corforward = cleanforward.replace(".clean.fq", ".clean.cor.fq")
    if not os.path.isfile(corforward):
        # Define /dev/null
        FNULL = open(os.devnull, 'wb')
        if not os.path.isfile(countsfile):
            # Define the command
            if cleanreverse:
                jellycommand = "jellyfish count -C -m 15 -s 10M -t 12 -o %s %s %s" \
                               % (countsfile, cleanforward, cleanreverse)
            else:
                jellycommand = "jellyfish count -C -m 15 -s 10M -t 12 -o %s %s" % (countsfile, cleanforward)
            # Run the command
            subprocess.call(jellycommand, shell=True, stdout=FNULL, stderr=FNULL)
        if not os.path.isfile("%s.txt" % countsfile):
            dumpcommand = "jellyfish dump -c -o %s.txt %s" % (countsfile, countsfile)
            subprocess.call(dumpcommand, shell=True, stdout=FNULL, stderr=FNULL)
        pathname = os.path.split(countsfile)[0]
        if not os.path.isfile("%s/cutoff.txt" % pathname):
            cutoffcall = "cd %s && cov_model.py --int --path=%s %s.txt > %s_cutoffDump.txt" % (pathname, pathname, countsfile, pathname)
            subprocess.call(cutoffcall, shell=True, stdout=FNULL, stderr=FNULL)
        if os.path.isfile('%s/cutoff.txt' % pathname):
            # Open the cutoff file and read the first line
            cutoff = open('%s/cutoff.txt' % pathname).readline().rstrip()
            if cutoff == "Inf":
                cutoff = 1
        else:
            # Otherwise set the cutoff to 1 if the cutoff command failed
            cutoff = 1
        quakeCorrect = "correct -f %s/fastqFiles.txt -k 15 -m %s/counts.txt -c %s -p 24" % (pathname, pathname, cutoff)
        # Run the correct command
        subprocess.call(quakeCorrect, shell=True, stdout=FNULL, stderr=FNULL)
    # Move the trimmed and corrected files to the appropriate folder
    if correverse:
        for fastqfile in [corforward, correverse]:
            finalname = os.path.basename(fastqfile).replace(".clean.cor.fq", ".fastq")
            # As the geneSippr will move these files into folders once it begins processing, get the folder names,
            # so the script won't try to copy over the files each time
            foldername = finalname.split("_R")[0]
            # Only try to copy the files if they have not already been copied over
            if not os.path.isfile("%s/%s" % (seqFolder, finalname)) and not os.path.isdir(foldername):
                shutil.copyfile(fastqfile, "%s/%s" % (seqFolder, finalname))
    else:
        finalname = os.path.basename(corforward).replace(".clean.cor.fq", ".fastq")
        # As the geneSippr will move these files into folders once it begins processing, get the name of the folders,
        # so the script won't try to copy over the files each time
        foldername = finalname.split("_R")[0]
        # Only try to copy the files if they have not already been copied over
        if not os.path.isfile("%s/%s" % (seqFolder, finalname)) and not os.path.isdir(foldername):
            shutil.copyfile(corforward, "%s/%s" % (seqFolder, finalname))

# Add an argument parser if this script is called as a stand-alone
if __name__ == "__main__":
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser

    #Parser for arguments
    parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s v1.0')
    parser.add_argument('-p', '--path', required=True, help='Specify input directory')
    parser.add_argument('-s', '--sequencePath', required=True, help='Path of .fastq(.gz) files to process. If not '
                        'provided, the default path of "path/sequences" will be used')
    parser.add_argument('-t', '--targetPath', required=False, help='Path of target files to process. If not '
                        'provided, the default path of "path/targets" will be used')
    parser.add_argument('-m', '--miSeqPath', required=False, help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', '--miSeqFolder', required=False, help='Name of the folder containing MiSeq run data')
    parser.add_argument('-r1', '--readLengthForward', required=False, help='Length of forward reads to use. Can specify'
                        '"full" to take the full length of forward reads specified on the SampleSheet')
    parser.add_argument('-r2', '--readLengthReverse', required=False, default=0, help='Length of reverse reads to use. '
                        'Can specify "full" to take the full length of reverse reads specified on the SampleSheet')
    parser.add_argument('-c', '--customSampleSheet', required=False, help='Path of folder containing a custom sample '
                        'sheet (still must be named "SampleSheet.csv")')
    parser.add_argument('-P', '--projectName', required=False, help='A name for the analyses. If nothing is provided, then '
                        'the "Sample_Project" field in the provided sample sheet will be used. Please note that bcl2fastq '
                        'creates subfolders using the project name, so if multiple names are provided, the results will be '
                        'split as into multiple projects')


    # Get the arguments into a list
    args = vars(parser.parse_args())

    # Define variables from the arguments - there may be a more streamlined way to do this
    path = os.path.join(args['path'], "")
    miseqPath = os.path.join(args['miSeqPath'], "")
    customSampleSheet = os.path.join(args['customSampleSheet'], "")

    # If these variables are not defined, set to default values
    if args['sequencePath']:
        sequencePath = os.path.join(args['sequencePath'], "")
    else:
        sequencePath = "%ssequences/" % path

    if args['targetPath']:
        targetPath = os.path.join(args['targetPath'], "")
    else:
        targetPath = "%stargets/" % path

    targets = glob("%s*.fasta" % targetPath)
    miSeqFolder = args['miSeqFolder']
    projectName = args['projectName']
    forwardreads = args['readLengthForward']
    reversereads = args['readLengthReverse']

    # Run the bcl2fastq process
    createFastq(miseqPath, miSeqFolder, path, projectName, forwardreads, reversereads, customSampleSheet)

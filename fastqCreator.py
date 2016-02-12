import errno
import os
import re
import shutil
import subprocess
import time
from glob import glob
from multiprocessing import Pool

__author__ = 'adamkoziol'


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    try:
        os.makedirs(inpath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def createfastq(miseqpath, miseqfolder, path, project, forwardreads, reversereads, customsamplesheet):
    """finds the most recent run folder on the MiSeq and determines which cycle the run is currently on. Runs bcl2fastq
    to create .fastq files with the number of reads specified in SampleSheet.csv
    :param customsamplesheet: (optional) name and path of custom sample sheet
    :param reversereads: (optional) length of forwards reads to use
    :param forwardreads: (optional) length of reverse reads to use
    :param project: (optional) name of the project
    :param path: string of the path of the analyses
    :param miseqfolder: name of the folder to use in the analyses
    :param miseqpath: path of the MiSeq analyses
    """
    if not miseqfolder:
        # Get the folders into a list
        miseqfolders = glob("%s/*/" % miseqpath)
        # The folder on interest is the most recent one, so sort the folders and choose the last on in the list
        miseqfolder = sorted(miseqfolders)[-1]
        foldername = miseqfolder.split("/")[-2]
    else:
        # Add the path to the folder name
        foldername = miseqfolder
        miseqfolder = miseqpath + miseqfolder

    # Get the flow cell id from the folder name - it should be the last part of the name, and will be separated by a "-"
    try:
        flowcellid = miseqfolder.split("_")[1]
    except IndexError:
        flowcellid = "NA"
    # Initialise the folder in which to place the .fastq sequences once they are created
    seqfolder = path + "sequences"
    folderpath = seqfolder + "/" + foldername
    # Make the folderPath path if necessary
    make_path(folderpath)
    # Initialise the reads list to be used for storing the length of the forward and reverse reads
    reads = []
    # Initialise the variables to store data from the original sample sheet
    investigator = ""
    forwardlength = ""
    reverselength = ""
    projectname = ""
    # bcl2fastq requires an older version of the sample sheet, this recreates the required version
    # Create the new sample sheet
    with open("%s/%s/SampleSheet_modified.csv" % (seqfolder, foldername), "wb") as modifiedSampleSheet:
        # Write the required headings to the file
        modifiedSampleSheet.write(
            "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n")
        # Open the current sample sheet - if there is a custom sample sheet provided, use it
        if customsamplesheet:
            originalsamplesheet = open("%sSampleSheet.csv" % customsamplesheet, "rb")
        else:
            originalsamplesheet = open("%s/SampleSheet.csv" % miseqfolder, "rb")
        for line in originalsamplesheet:
            # As this is a csv file, data are separated by commas - split on commas
            data = line.split(",")
            if "Investigator Name" in line:
                investigator = data[1].rstrip()
            # Here's a Perl-like solution for reading lines after a regex match
            # Perform the regex
            elif "[Reads]" in line:
                # Now read all sublines going forward
                for subline in originalsamplesheet:
                    # Stop reading once "Settings" is matched
                    if "[Settings]" in subline:
                        break
                    reads.append(subline.rstrip())
                # Grab the number of reads in the forward and reverse direction
                forwardlength = int(reads[0].replace(",", "").rstrip())
                reverselength = int(reads[1].replace(",", "").rstrip())
            # Get the data for each sample from the sheet
            elif "Sample_ID" in line:
                # Iterate through the lines following "Sample_ID"
                for subline in originalsamplesheet:
                    # Split on the ,
                    subdata = subline.split(",")
                    # Populate
                    sampleid = subdata[0]
                    index = "%s-%s" % (subdata[5], subdata[7])
                    # Not fully tested - attempts to catch sample sheets without supplied Sample_Project and/or
                    # Descriptions
                    # Only run this if the project variable has not been provided
                    if not project:
                        # If the there is a project name in the sample sheet
                        if subdata[8]:
                            projectname = subdata[8].rstrip()
                        # Otherwise supply a NA
                        else:
                            projectname = "NA"
                    else:
                        projectname = project
                    # Same as above look for the data, and provide a NA if it's not provided
                    if subdata[9]:
                        description = subdata[9].rstrip()
                    else:
                        description = "NA"
                    # Write the data to the modified sample sheet
                    modifiedSampleSheet.write("%s,1,%s,no_ref,%s,%s,N,NA,%s,%s\n" % (flowcellid, sampleid, index,
                                                                                     description, investigator,
                                                                                     projectname))
    # Calculate the number of cycles completed
    cycles = glob("%s/Data/Intensities/BaseCalls/L001/C*" % miseqfolder)
    # As the length of reads can be specified as "full", use the forwardlength and reverselength variables to populate
    # the forwardreads and reversereads as required
    if forwardreads == "full":
        forwardreads = forwardlength
    else:
        forwardreads = int(forwardreads)
    if reversereads == "full":
        reversereads = reverselength
    else:
        reversereads = int(reversereads)
    # As the number of cycles required is the number of forward reads + the index(8) + the second index(8)
    # Also set the basemask variable as required
    if forwardreads:
        if reversereads:
            readsneeded = forwardreads + reversereads
            basemask = "Y%sn*,I8,I8,Y%sn*" % (forwardreads, reversereads)
            nohup = "nohup make -j 16"
        else:
            readsneeded = forwardreads + 16
            basemask = "Y%sn*,I8,I8,n*" % forwardreads
            nohup = "nohup make -j 16 r1"
    else:
        if reverselength:
            readsneeded = forwardlength + reverselength
            basemask = "Y%sn*,I8,I8,Y%sn*" % (forwardlength, reverselength)
            nohup = "nohup make -j 16"
        else:
            readsneeded = forwardlength + 16
            basemask = "Y%sn*,I8,I8,n*" % forwardlength
            nohup = "nohup make -j 16 r1"

    # A while loop to calculate when there are enough cycles to perform the analysis as specified in the arguments
    while len(cycles) < readsneeded:
        print "Currently at %s cycles. Need to wait until the MiSeq reaches cycles %s" % (len(cycles), readsneeded)
        time.sleep(30)
        cycles = glob("%s/Data/Intensities/BaseCalls/L001/C*" % miseqfolder)
    if not os.path.isdir("%s/Project_%s" % (folderpath, projectname)):
        # Call configureBclToFastq.pl
        print "Running bcl2fastq"
        bclcall = "configureBclToFastq.pl --input-dir %s/Data/Intensities/BaseCalls/ " \
                  "--output-dir %s --force --sample-sheet %s/SampleSheet_modified.csv " \
                  "--mismatches 1 --no-eamss --fastq-cluster-count 0 --compression none --use-bases-mask %s" \
                  % (miseqfolder, folderpath, folderpath, basemask)
        print bclcall
        # Define /dev/null
        fnull = open(os.devnull, 'wb')
        # Run the commands
        subprocess.call(bclcall, shell=True, stdout=fnull, stderr=fnull)
        nohupcall = "cd %s/%s && %s" % (seqfolder, foldername, nohup)
        subprocess.call(nohupcall, shell=True, stdout=fnull, stderr=fnull)

    # Find all the files created using a nested list comprehension
    # Firstly, list all the directories using os.listdir (the second part of the comprehension)
    # Then look for .gz files in each of the newly discovered folders
    gzfiles = [[sorted(glob("%s/Project_%s/%s/*.gz" % (folderpath, projectname, f)))]
               for f in sorted(os.listdir("%s/Project_%s" % (folderpath, projectname)))]

    # Initialise lists to hold arguments for multiprocessing modules
    bbdukargs = []
    quakeargs = []
    # Iterate through the nested lists
    for i in range(len(gzfiles)):
        # Second level
        for j in range(len(gzfiles[i])):
            # Third level
            for k in range((len(gzfiles[i][j]))):
                # File names are strings of the kth element of the jth element of the ith element of gzfiles
                filename = str(gzfiles[i][j][k])
                # Substitute in the sequence number for the indexes used. Should look something like:
                # ..../Sample_2015-SEQ-1106/2015-SEQ-1106_GTAGAGGA-CTAAGCCT_L001_R1_001.fastq.gz becomes
                # Sample_2015-SEQ-1106/2015-SEQ-1106_S1_L001_R1_001.fastq.gz
                subname = re.sub("\w{8}-\w{8}", "S%s" % str(i + 1), filename)
                # Rename the files appropriately
                os.rename(filename, subname)
            # Treat paired vs. unpaired reads differently
            forward = ""
            reverse = ""
            cleanreverse = ""
            if len(gzfiles[i][j]) == 2:
                forward, reverse = gzfiles[i][j]
                # Uses lambda to substitute out the "_S1_L001_" and replace the "_001.fastq.gz" portions of the filename
                cleanforward, cleanreverse = map(lambda x: re.sub("_S\d+_L001", "",
                                                                  x.replace("_001.fastq.gz", ".clean.fq")),
                                                 gzfiles[i][j])
            elif len(gzfiles[i][j]) == 1:
                forward = gzfiles[i][j][0]
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
        # fnull = open(os.devnull, 'wb')
        # Run the command - using os.system right now, as subprocess.call failed the last time I tried it
        os.system(bbdukcall)
        # subprocess.call(bbdukcall, shell=True, stdout=fnull, stderr=fnull)


def quake((cleanforward, cleanreverse, countsfile, seqfolder)):
    """Runs jellyfish count with the supplied variables in a multiprocessed fashion"""
    correverse = ""
    if cleanreverse:
        # Define the names of the corrected fastq files
        corforward, correverse = map(lambda x: x.replace(".clean.fq", ".clean.cor.fq"), [cleanforward, cleanreverse])
    else:
        corforward = cleanforward.replace(".clean.fq", ".clean.cor.fq")
    if not os.path.isfile(corforward):
        # Define /dev/null
        fnull = open(os.devnull, 'wb')
        if not os.path.isfile(countsfile):
            # Define the command
            if cleanreverse:
                jellycommand = "jellyfish count -C -m 15 -s 10M -t 12 -o %s %s %s" \
                               % (countsfile, cleanforward, cleanreverse)
            else:
                jellycommand = "jellyfish count -C -m 15 -s 10M -t 12 -o %s %s" % (countsfile, cleanforward)
            # Run the command
            subprocess.call(jellycommand, shell=True, stdout=fnull, stderr=fnull)
        if not os.path.isfile("%s.txt" % countsfile):
            dumpcommand = "jellyfish dump -c -o %s.txt %s" % (countsfile, countsfile)
            subprocess.call(dumpcommand, shell=True, stdout=fnull, stderr=fnull)
        pathname = os.path.split(countsfile)[0]
        if not os.path.isfile("%s/cutoff.txt" % pathname):
            cutoffcall = "cd %s && cov_model.py --int --path=%s %s.txt > %s_cutoffDump.txt" % (
                pathname, pathname, countsfile, pathname)
            subprocess.call(cutoffcall, shell=True, stdout=fnull, stderr=fnull)
        if os.path.isfile('%s/cutoff.txt' % pathname):
            # Open the cutoff file and read the first line
            cutoff = open('%s/cutoff.txt' % pathname).readline().rstrip()
            if cutoff == "Inf":
                cutoff = 1
        else:
            # Otherwise set the cutoff to 1 if the cutoff command failed
            cutoff = 1
        quakecorrect = "correct -f %s/fastqFiles.txt -k 15 -m %s/counts.txt -c %s -p 24" % (pathname, pathname, cutoff)
        # Run the correct command
        subprocess.call(quakecorrect, shell=True, stdout=fnull, stderr=fnull)
    # Move the trimmed and corrected files to the appropriate folder
    if correverse:
        for fastqfile in [corforward, correverse]:
            finalname = os.path.basename(fastqfile).replace(".clean.cor.fq", ".fastq")
            # As the geneSippr will move these files into folders once it begins processing, get the folder names,
            # so the script won't try to copy over the files each time
            foldername = finalname.split("_R")[0]
            # Only try to copy the files if they have not already been copied over
            if not os.path.isfile("%s/%s" % (seqfolder, finalname)) and not os.path.isdir(foldername):
                shutil.copyfile(fastqfile, "%s/%s" % (seqfolder, finalname))
    else:
        finalname = os.path.basename(corforward).replace(".clean.cor.fq", ".fastq")
        # As the geneSippr will move these files into folders once it begins processing, get the name of the folders,
        # so the script won't try to copy over the files each time
        foldername = finalname.split("_R")[0]
        # Only try to copy the files if they have not already been copied over
        if not os.path.isfile("%s/%s" % (seqfolder, finalname)) and not os.path.isdir(foldername):
            shutil.copyfile(corforward, "%s/%s" % (seqfolder, finalname))


# Add an argument parser if this script is called as a stand-alone
if __name__ == "__main__":
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser

    # Parser for arguments
    parser = ArgumentParser(description='Perform modelling of parameters for GeneSipping')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s v1.0')
    parser.add_argument('-p', '--path', required=True, help='Specify input directory')
    parser.add_argument('-s', '--sequencePath', required=False, help='Path of .fastq(.gz) files to process. If not '
                        'provided, the default path of "path/sequences" will be used')
    parser.add_argument('-t', '--targetPath', required=False, help='Path of target files to process. If not '
                        'provided, the default path of "path/targets" will be used')
    parser.add_argument('-m', '--miseqpath', required=False, help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', '--miseqfolder', required=False, help='Name of the folder containing MiSeq run data')
    parser.add_argument('-r1', '--readLengthForward', required=False, help='Length of forward reads to use. Can specify'
                        '"full" to take the full length of forward reads specified on the SampleSheet')
    parser.add_argument('-r2', '--readLengthReverse', required=False, default=0, help='Length of reverse reads to use. '
                        'Can specify "full" to take the full length of reverse reads specified on the SampleSheet')
    parser.add_argument('-c', '--customsamplesheet', required=False, help='Path of folder containing a custom sample '
                        'sheet (still must be named "SampleSheet.csv")')
    parser.add_argument('-P', '--projectname', required=False, help='A name for the analyses. If nothing is provided, '
                        'then the "Sample_Project" field in the provided sample sheet will be used. Please note that '
                        'bcl2fastq creates subfolders using the project name, so if multiple names are provided, the '
                        'results will be split into multiple projects')

    # Get the arguments into a list
    args = vars(parser.parse_args())
    # Define variables from the arguments - there may be a more streamlined way to do this
    patharg = os.path.join(args['path'], "")
    miseqPath = os.path.join(args['miseqpath'], "")
    if args['customsamplesheet']:
        customsamplesheetarg = os.path.join(args['customsamplesheet'], "")
    else:
        customsamplesheetarg = ""
    # If these variables are not defined, set to default values
    if args['sequencePath']:
        sequencePath = os.path.join(args['sequencePath'], "")
    else:
        sequencePath = "%ssequences/" % patharg
    if args['targetPath']:
        targetPath = os.path.join(args['targetPath'], "")
    else:
        targetPath = "%stargets/" % patharg

    targets = glob("%s*.fasta" % targetPath)
    miseqfolderarg = args['miseqfolder']
    projectnamearg = args['projectname']
    forwardreadsarg = args['readLengthForward']
    reversereadsarg = args['readLengthReverse']

    # Run the bcl2fastq process
    createfastq(miseqPath, miseqfolderarg, patharg, projectnamearg, forwardreadsarg, reversereadsarg,
                customsamplesheetarg)

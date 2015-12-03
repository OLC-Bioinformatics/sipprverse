import json, os, sys, errno, time, re, operator
from collections import defaultdict


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


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)


def alleleSplitter(alleleNames):
    # Multiple try-excepts. Maybe overly complicated, but I couldn't get it work any other way
    # This (hopefully) accounts for all the possible naming schemes for the alleles
    try:  # no split - just allele numbers e.g. >12
        match = re.search(r"(>\d+)", alleleNames)
        alleleNumber = match.group().split(">")[1]
        allelePreNumber = ""
    except (IndexError, AttributeError):
        try:  # split on "_" e.g. >AROC_12
            # alleleNumber is the number of the allele(!). It should be different for each allele
            alleleNumber = alleleNames.split("_")[1]
            # allelePreNumber is anything before the allele number. It should be the same for each allele
            allelePreNumber = alleleNames.split("_")[0]
        except IndexError:
            try:  # split on "-" e.g. >AROC-12
                alleleNumber = alleleNames.split("-")[1]
                allelePreNumber = alleleNames.split("-")[0]
            except IndexError:
                try:  # split on " " e.g. >AROC 12
                    alleleNumber = alleleNames.split(" ")[1]
                    allelePreNumber = alleleNames.split(" ")[0]
                except IndexError:
                    try:  # split on change from letters to numbers e.g. >AROC12
                        match = re.search(r"(>[A-z/a-z]+)(\d+)", alleleNames)
                        alleleNumber = match.groups()[1]
                        allelePreNumber = match.groups()[0]
                    except (IndexError, AttributeError):
                        alleleNumber = alleleNames
                        allelePreNumber = alleleNames
    # Return the variables
    return int(alleleNumber), allelePreNumber


def profilR(profileFile, profileType):
    """Creates a dictionary from the profile scheme"""
    # Initialise the dictionary
    profileData = defaultdict(make_dict)
    geneDict = {}
    profileDict = {}
    for strain in profileFile:
        lastEntry = ""
        geneList = []
        geneSet = set()
        if profileFile[strain][profileType]:
            # print strain, profileFile[strain]["MLSTprofile"]
            # The gene names are present in the first line of the profile file
            # Note: if the profile files are ever updated, then the clonal complex column must be removed
            # Make the json filename from profileFile - it might already be .json, but this doesn't take long to do
            JSONProfile = "%s.json" % os.path.splitext(profileFile[strain][profileType])[0]
            # If this scheme has previously been used, then the profileData dictionary is written to disk for increased speed.
            # Parsing a json file was approximately 10 times faster than parsing the original tab-delimited file
            # Get the MLST profiles for each sequence type
            # Don't do this if the .json profile has previously been created
            if not os.path.isfile(JSONProfile):
                with open(profileFile[strain][profileType]) as profile:
                    # Files have to be in tab-delimited format
                    header = profile.readline().rstrip().split("\t")
                    # As certain typing schemes are not in any discernible order, using a naturally ordered list instead of a
                    # dictionary to store the names is a good idea
                    for gene in header:
                        # The first column must have "ST" in the header
                        if not "ST" in gene:
                            dotter()
                            geneList.append(gene)
                    for line in profile:
                        # Grab the name of the last profile
                        # MLSTcount will used to associate the gene name in header to the allele (e.g. adk 12)
                        MLSTcount = 1
                        # Don't need to get the header information again
                        # if not "ST" in line:
                            # len(header) will be the same length as the data in line
                        while MLSTcount < len(header):
                            # Remove newlines and split on tabs
                            data = line.rstrip().split("\t")
                            # Populate profileData with the sequence type, gene name, and the allele number
                            profileData[data[0]][header[MLSTcount]] = data[MLSTcount]
                            # Increment
                            MLSTcount += 1
                    # Split the name (if necessary) to just have the profile number
                    # Write the json file to disk
                    JSONreport = open(JSONProfile, "wb")
                    output = json.dumps(profileData, sort_keys=True, indent=4, separators=(',', ': '))
                    JSONreport.write(output)
                    JSONreport.close()
            else:
                with open(JSONProfile, "rb") as jsonReport:
                    # Load the data
                    profileData = json.load(jsonReport)
                    for profile in profileData:
                        for gene in profileData[profile]:
                            geneSet.add(gene)
                geneList = list(geneSet)
            geneDict[strain] = geneList
            profileDict[strain] = profileData
        else:
            geneDict[strain] = "N"
            profileDict[strain] = "N"
        # print "\n".join(geneList)
    # print json.dumps(geneDict, sort_keys=True, indent=4, separators=(',', ': '))
    dotter()
    return profileDict, geneDict


def sequenceTyper(MLSTmatches, profileDict, geneDict):
    """Determines the sequence type of each strain based on comparisons to sequence type profiles"""
    # global MLSTseqType
    # global bestDict
    # Initialise variables
    bestDict = defaultdict(make_dict)
    MLSTseqType = defaultdict(make_dict)
    resultProfile = defaultdict(make_dict)
    for genome in profileDict:
        if MLSTmatches[genome]:
            # print genome, "\t".join(MLSTmatches[genome].keys())
            header = 0
            alleleCount = 0
            multiAllele = []
            multiPercent = []
            bestMatch = defaultdict(int)
            bestCount = 0
            # Iterate through the genomes
            # for genome in plusdict:
            # global resultProfile
            # Initialise bestMatch[genome] with an integer - this will eventually be replaced by the number of matches
            bestMatch[genome] = defaultdict(int)
            # For each gene in plusdict[genome]
            for gene in MLSTmatches[genome]:
                # print genome, gene
                # Clear the appropriate count and lists
                alleleCount = 0
                multiAllele = []
                multiPercent = []
                multiFold = []
                for allele in MLSTmatches[genome][gene]:
                    for percentID, foldCoverage in MLSTmatches[genome][gene][allele].iteritems():
                        # print genome, gene, allele, percentID, foldCoverage
                        # "N" alleles screw up the allele splitter function
                        if allele != "N":
                            # Use the alleleSplitter function to get the allele number
                            alleleNumber, allelePreNumber = alleleSplitter(allele)
                            # print genome, gene, allele, allelePreNumber, alleleNumber, percentID, foldCoverage
                            # Append as appropriate - alleleNumber is treated as an integer for proper sorting
                            multiAllele.append(int(alleleNumber))
                            multiPercent.append(percentID)
                            multiFold.append(foldCoverage)
                        # If the allele is "N"
                        else:
                            # Append "N" and a percent identity of 0
                            multiAllele.append("N")
                            multiPercent.append(0)
                            multiFold.append(0)
                    # Trying to catch cases that where the allele isn't "N", but can't be parsed by alleleSplitter
                    if not multiAllele:
                        multiAllele.append("N")
                        multiPercent.append(0)
                        multiFold.append(0)
                    # Populate bestDict with genome, gene, alleles - joined with a space (this was written like this because
                    # allele is a list generated by the .iteritems() above, the percent identity, and the fold coverage
                    bestDict[genome][gene][" ".join(str(allele) for allele in sorted(multiAllele))][multiPercent[0]] = multiFold[0]
                # Find the profile with the most alleles in common with the query genome
                for sequenceType in profileDict[genome]:
                    # Reset counts to 0
                    matchCount = 0
                    bestCount = 0
                    # The number of genes in the analysis
                    header = len(profileDict[genome][sequenceType])
                    # refallele is the allele number of the sequence type
                    refAllele = profileDict[genome][sequenceType][gene]
                    # If there are multiple allele matches for a gene in the reference profile e.g. 10 692
                    # if len(refAllele.split(" ")) > 1:
                    #     # Map the split (on a space) alleles as integers - if they are treated as integers,
                    #     # the alleles will sort properly
                    #     intRefAllele = map(int, refAllele.split(" "))
                    #     # Create a string of the joined, sorted allleles
                    #     sortedRefAllele = " ".join(str(allele) for allele in sorted(intRefAllele))
                    # else:
                    #     # Use the reference allele as the sortedRefAllele
                    #     sortedRefAllele = refAllele
                    #     # print genome, sequenceType, refAllele
                    refAlleleList = refAllele.split(" ")
                    for allele in bestDict[genome][gene]:
                        # If the allele for the gene in the query genome matches the allele in the reference profile,
                        # add the result to the bestMatch dictionary. Because genes with multiple alleles were sorted the
                        # same way, these strings with multiple alleles will match e.g. 10 692 will never be 692 10
                        for rAllele in refAlleleList:
                            if allele == rAllele:
                                # Increment the number of matches to each profile
                                # print genome, gene, allele, rAllele, sortedRefAlleleList
                                bestMatch[genome][sequenceType] += 1
                                break
            # Get the best number of matches
            # From: https://stackoverflow.com/questions/613183/sort-a-python-dictionary-by-value
            sortedMatches = sorted(bestMatch[genome].items(), key=operator.itemgetter(1), reverse=True)[0][1]
            # If there are fewer matches than the total number of genes in the typing scheme
            if int(sortedMatches) < header:
                # Iterate through the sequence types and the number of matches in bestDict for each genome
                for sequenceType, matches in bestMatch[genome].iteritems():
                    # If the number of matches for a profile matches the best number of matches
                    if matches == sortedMatches:
                        # Iterate through the gene in the analysis
                        for gene in profileDict[genome][sequenceType]:
                            # Get the reference allele as above
                            refAllele = profileDict[genome][sequenceType][gene]
                            # As above get the reference allele split and ordered as necessary
                            if len(refAllele.split(" ")) > 1:
                                intRefAllele = map(int, refAllele.split(" "))
                                sortedRefAllele = " ".join(str(allele) for allele in sorted(intRefAllele))
                            else:
                                sortedRefAllele = refAllele
                            # If there are missing alleles, populate the dictionary to indicate this
                            refAlleleList = refAllele.split(" ")
                            try:
                                # Somewhat complicated here - dictionary is genome, sequence type,
                                # number of matches to sequence type, gene, observed allele - taken from bestDict,
                                # reference allele, percent identity taken from bestDict - the keys of the values of bestDict[genome][gene]
                                # fold coverage - taken from bestDict
                                # MLSTseqType[genome][sequenceType][sortedMatches][gene][str(bestDict[genome][gene].keys()[0])][sortedRefAllele][str(bestDict[genome][gene].values()[0].keys()[0])] = str(bestDict[genome][gene].values()[0].values()[0])
                                for rAllele in refAlleleList:
                                    # if allele == rAllele:
                                    if rAllele == str(bestDict[genome][gene].keys()[0]):
                                        MLSTseqType[genome][sequenceType][sortedMatches][gene][str(bestDict[genome][gene].keys()[0])][str(bestDict[genome][gene].values()[0].keys()[0])] = str(bestDict[genome][gene].values()[0].values()[0])
                                        break
                                if not MLSTseqType[genome][sequenceType][sortedMatches][gene][str(bestDict[genome][gene].keys()[0])][str(bestDict[genome][gene].values()[0].keys()[0])]:
                                    MLSTseqType[genome][sequenceType][sortedMatches][gene].clear()
                                    MLSTseqType[genome][sequenceType][sortedMatches][gene]["%s (%s)" % (str(bestDict[genome][gene].keys()[0]), sortedRefAllele)][str(bestDict[genome][gene].values()[0].keys()[0])] = str(bestDict[genome][gene].values()[0].values()[0])
                            except IndexError:
                                MLSTseqType[genome][sequenceType][sortedMatches][gene]["N (%s)" % sortedRefAllele][0] = 0
            #     # Add the new profile to the profile file (if the option is enabled)
            #     if updateProfile:
            #         reProfilR(int(header), profileFile, geneList, genome)
            # Otherwise, the query profile matches the reference profile
            else:
                # Iterate through best match
                for sequenceType, matches in bestMatch[genome].iteritems():
                    # print genome, sequenceType
                    if matches == sortedMatches:
                        for gene in profileDict[genome][sequenceType]:
                            for allele in MLSTmatches[genome][gene]:
                                # Use the alleleSplitter function to get the allele number
                                alleleNumber, allelePreNumber = alleleSplitter(allele)
                                for percentID, foldCoverage in MLSTmatches[genome][gene][str(allele)].iteritems():
                                    # Populate resultProfile with the genome, best match to profile, number of
                                    # matches to the profile, gene, query allele(s), reference allele(s), percent identity, and fold coverage
                                    resultProfile[genome][sequenceType][sortedMatches][gene][alleleNumber][percentID] = foldCoverage
        dotter()
    resultProfile.update(MLSTseqType)
    # print json.dumps(resultProfile, sort_keys=True, indent=4, separators=(',', ': '))
    return resultProfile


def MLSTreportMaker(seqDict, sequenceTypes, analysisType, reportFolder, organismDict, organismList, path):
    """Helpful comment"""
    completeheader = ""
    completestring = ""
    for organism in organismList:
        compiledResultString = ""
        csvheader = ""
        orgcount = 0
        for strain in organismDict:
            # reportFolder = "%sreports/%s" % (path, time.strftime("%Y.%m.%d.%H.%M.%S"))
            # make_path(reportFolder)
            # compiledResultString = ""
            # csvheader = ""
            for currentorganism in organismDict[strain]:
                if organism == currentorganism:
                    orgcount += 1
                    # Create variables as required
                    baitType = seqDict[strain]["bait"]["fastqFiles"].keys()[0]
                    fastqDir = os.path.split(seqDict[strain]["bait"]["fastqFiles"][baitType][0])[0]
                    reportDir = "%s/reports" % fastqDir
                    reportName = "%s/%s_%s_reports.csv" % (reportDir, strain, baitType)
                    make_path(reportDir)
                    csvfile = open(reportName, "wb")
                    # Get the header started
                    csvheader = "Strain,SequenceType,Matches,"

                    # Initialise the headerGenes variable
                    headerGenes = ''
                    for target in sorted(seqDict[strain]["targets"][analysisType]):
                        targetName = os.path.basename(target).split(".")[0]
                        # print strain, reportDir, targetName
                        # Append each gene to headerGenes
                        headerGenes += "%s," % targetName
                    # Append headerGenes to the header
                    csvheader += headerGenes
                    csvheader += "\n"
                    if orgcount == 1:
                        completestring += "Strain,Genus,SequenceType,Matches,"
                        completestring += headerGenes
                        completestring += "\n"
                    # Write the header to the report
                    csvfile.write(csvheader)
                    # Reset resultString to a newline
                    resultString = ""
                    # Reset the count to 0
                    resProfileCount = 0
                    # Iterate through the sequence types
                    for sequenceType in sequenceTypes[strain]:
                        if resProfileCount == 0:
                            # Put the genome and the sequence type in the result string
                            resultString += "%s,%s," % (strain, sequenceType)
                            compiledResultString += "%s,%s," % (strain, sequenceType)
                            # if orgcount == 1:
                            completestring += "%s,%s,%s," % (strain, organism, sequenceType)
                            # else:
                            #     completestring += "%s,,%s," % (strain, sequenceType)
                        else:
                            resultString += ",%s," % sequenceType
                            compiledResultString += ",%s," % sequenceType
                            # if orgcount == 1:
                            completestring += ",,%s," % sequenceType
                            # else:
                            #     completestring += ",,,%s," % sequenceType

                        # Report the number of matches to the profile
                        for numMatches in sequenceTypes[strain][sequenceType]:
                            # Append these matches
                            resultString += "%s," % numMatches
                            compiledResultString += "%s," % numMatches
                            completestring += "%s," % numMatches
                            # Go through the ordered list of genes
                            for target in seqDict[strain]["targets"][analysisType]:
                                targetName = os.path.basename(target).split(".")[0]
                                # Add each allele to the result string
                                for allele in sequenceTypes[strain][sequenceType][numMatches][targetName]:
                                    # Add the allele to resultString
                                    resultString += "%s," % allele
                                    compiledResultString += "%s," % allele
                                    completestring += "%s," % allele
                                    # Increment the count variable
                                    resProfileCount += 1
                        resultString += "\n"
                        compiledResultString += "\n"
                        completestring += "\n"
                    # Append the resultString to the report
                    csvfile.write(resultString)
                    csvfile.close()


            compiledcsvfile = open("%s/%s_%s_results.csv" % (reportFolder, analysisType, organism), "wb")
            compiledcsvfile.write(csvheader)
            compiledcsvfile.write(compiledResultString)
            compiledcsvfile.close()
    completecsvfile = open("%s/%s_complete_results.csv" % (reportFolder, analysisType), "wb")
    # completecsvfile.write(completeheader)
    completecsvfile.write(completestring)
    completecsvfile.close()


def rMLSTreportMaker(seqDict, sequenceTypes, analysisType, reportFolder, path):
    """Helpful comment"""
    completeheader = ""
    completestring = ""
    compiledResultString = ""
    csvheader = ""
    header = True
    for strain in seqDict:
        # reportFolder = "%sreports/%s" % (path, time.strftime("%Y.%m.%d.%H.%M.%S"))
        # make_path(reportFolder)
        # csvheader = ""
        # Create variables as required
        baitType = seqDict[strain]["bait"]["fastqFiles"].keys()[0]
        fastqDir = os.path.split(seqDict[strain]["bait"]["fastqFiles"][baitType][0])[0]
        reportDir = "%s/reports" % fastqDir
        reportName = "%s/%s_%s_reports.csv" % (reportDir, strain, baitType)
        make_path(reportDir)
        csvfile = open(reportName, "wb")
        # Get the header started
        # Initialise the headerGenes variable
        headerGenes = ''
        if header:
            csvheader = "Strain,SequenceType,Matches,"

            for target in sorted(seqDict[strain]["targets"][analysisType]):
                targetName = os.path.basename(target).split(".")[0]
                # print strain, reportDir, targetName
                # Append each gene to headerGenes
                headerGenes += "%s," % targetName
            # Append headerGenes to the header
            csvheader += headerGenes
            # csvheader += "\n"
            completestring += "Strain,SequenceType,Matches,"
            completestring += headerGenes
            # completestring += "\n"
            # Write the header to the report
            csvfile.write(csvheader)
            header = False
        # Reset resultString to a newline
        resultString = ""
        # Reset the count to 0
        resProfileCount = 0
        # Iterate through the sequence types
        for sequenceType in sequenceTypes[strain]:
            if resProfileCount == 0:
                # Put the genome and the sequence type in the result string
                resultString += "\n%s,%s," % (strain, sequenceType)
                # if orgcount == 1:
                completestring += "\n%s,%s," % (strain, sequenceType)
                # else:
                #     completestring += "%s,,%s," % (strain, sequenceType)
            else:
                resultString += "\n,%s," % sequenceType
                # if orgcount == 1:
                completestring += "\n,%s," % sequenceType
                # else:
                #     completestring += ",,,%s," % sequenceType

            # Report the number of matches to the profile
            for numMatches in sequenceTypes[strain][sequenceType]:
                # Append these matches
                resultString += "%s," % numMatches
                completestring += "%s," % numMatches
                # Go through the ordered list of genes
                for target in seqDict[strain]["targets"][analysisType]:
                    targetName = os.path.basename(target).split(".")[0]
                    # Add each allele to the result string
                    for allele in sequenceTypes[strain][sequenceType][numMatches][targetName]:
                        # Add the allele to resultString
                        resultString += "%s," % allele
                        completestring += "%s," % allele
                        # Increment the count variable
                        resProfileCount += 1
        # Append the resultString to the report
        csvfile.write(resultString)
        csvfile.close()
    completestring += "\n"
    completecsvfile = open("%s/%s_complete_results.csv" % (reportFolder, analysisType), "wb")
    # completecsvfile.write(completeheader)
    completecsvfile.write(completestring)
    completecsvfile.close()


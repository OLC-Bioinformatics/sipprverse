import os
import json
import errno
import shutil
from collections import defaultdict


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    try:
        # os.makedirs makes parental folders as required
        os.makedirs(inpath)
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)


def commandr(samplenames, path):
    """
    Opens the *name*_metadataCollection.json file and extracts any commands performed from previous iterations of the
    pipeline on data in the supplied path
    :param samplenames: names of all the strains to process
    :param path: path of the folder containing the sequences and targets
    """
    # Initialise the command dictionary
    performedcommands = defaultdict(make_dict)
    # Open the *name*_metadataCollection.json file for each sample
    for name in samplenames:
        if os.path.isfile("%s/%s/%s_metadataCollection.json" % (path, name, name)):
            countsize = os.stat("%s/%s/%s_metadataCollection.json" % (path, name, name)).st_size
            if countsize != 0:
                with open("%s/%s/%s_metadataCollection.json" % (path, name, name)) as jsonreport:
                    # Load the data
                    jsondata = json.load(jsonreport)
                    if jsondata:
                        # Find the precise command used
                        for command in jsondata["commands"]:
                            # If the command exists, and is in the right format
                            if "N/A" not in jsondata["commands"][command] \
                                    and "defaultdict" not in jsondata["commands"][command]:
                                # Populate the performedCommands dictionary as appropriate
                                performedcommands[name][str(command)] = str(jsondata["commands"][command])
    # Return the dictionary
    return performedcommands


def jsonr(samplenames, path, metadata, filetype):
    """
    Creates a JSON report from a supplied metadata file
    :param samplenames: names of all the strains to process
    :param path: path of the folder containing the sequences and targets
    :param metadata: dictionary of metadata
    :param filetype: string of whether the json file is an intermediate "collections" or final file
    """
    # Make the reports folder as required
    reportpath = "%sreports" % path
    make_path(reportpath)
    for name in samplenames:
            # Create the .json file for each sample
            newpath = path + name
            reportname = "%s_metadata%s.json" % (name, filetype)
            jsonreport = open("%s/%s" % (newpath, reportname), "wb")
            # Print the JSON data to file
            output = json.dumps(metadata[name], sort_keys=True, indent=4, separators=(',', ': '))
            jsonreport.write(output)
            jsonreport.close()
            # Move all the reports to a common directory
            shutil.copy("%s/%s" % (newpath, reportname), "%s/%s" % (reportpath, reportname))

__author__ = 'mikeknowles'

from glob import glob
from re import match
from ARMICARD import decipher
import json, ARMISeekr




def blaster(path, targets, out, threshold, db, aro):
    if db == "both":
        db = ['ardb', 'card']
    else:
        db = [db]
    jsonfile = '%splusdict.json' % targets
    # if os.path.isfile(jsonfile):
    #     plusdict = json.load(open(jsonfile))
    #
    # else:
    markers = glob(path + "/*.fa*")
    for marker in markers:
        cardcheck = match("^\d{7}$", marker)
        if db == 'ardb' and cardcheck is not None:
            markers.remove(marker)
        elif db == 'card' and cardcheck is None:
            markers.remove(marker)

    plusdict = ARMISeekr.blaster(markers, targets, out, 'ARMI2')
    json.dump(plusdict, open(jsonfile, 'w'), sort_keys=True, indent=4, separators=(',', ': '))
    print json.dumps(plusdict, sort_keys=True, indent=4, separators=(',', ': '))
    antidict = json.load(open(aro))
    decipher(plusdict, antidict, out)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Antibiotic Resistance Marker Identifier:\n'
                                                 'Use to find markers for any bacterial genome')
    parser.add_argument('--version', action='version', version='%(prog)s v0.1')
    parser.add_argument('-i', '--input', required=True, help='Specify input fasta folder')
    parser.add_argument('-m', '--marker', required=True, help='Specify antibiotic markers folder')
    parser.add_argument('-o', '--output', required=True, help='Specify output folder for csv')
    parser.add_argument('-t', '--tab', type=str, required=True, help='tables file location')
    parser.add_argument('-c', '--cutoff', type=int, default=1, help='Threshold for maximum unique bacteria'
                                                                    ' for a single antibiotic')
    parser.add_argument('-d', '--db', default='both', help='Specify antibiotic markers database')
    args = vars(parser.parse_args())
    blaster(args['marker'], args['input'], args['output'], args['cutoff'], args['db'], args['tab'])
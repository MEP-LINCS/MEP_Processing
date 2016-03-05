import os
import sys

import synapseclient

# library(rGithubClient)
# synapseRawDataDir = "syn5706233"
synapseRawDataDir = "syn5706954"
synapseAnnotatedDataDir = "syn5706203"

def updateRow(row, key="name", path="./"):

    newrow = dict(row)

    filename = row[key]
    newrow.pop(key)

    return (os.path.abspath("%s/%s" % (path, filename)), newrow)

def process(row, syn, key, path, parentId):
    filename, annots = updateRow(row, key, path)
    f = synapseclient.File(filename, parentId=parentId)

    f = syn.store(f, forceVersion=False)

    syn.setAnnotations(f, annots)

def main():

    import csv
    import argparse

    parser = argparse.ArgumentParser("Upload raw data files.")
    parser.add_argument("-d", "--directory", type=str, help="Path to where raw data files are.")
    parser.add_argument("-k", "--key", type=str, default="name", help="Column in datafile with file name.")
    parser.add_argument("datafile", type=str, "csv file containing file name and annotations")

    args = parser.parse_args()

    print args

    syn = synapseclient.login(silent=True)

    with file(args.datafile) as f:
        reader = csv.DictReader(f)

        map(lambda x: process(x, syn, key=args.key, path=args.directory,
                              parentId=synapseRawDataDir),
            reader)

if __name__ == "__main__":
    main()

import os
import sys

import synapseclient

# synapseRawDataDir = "syn5706233"
# synapseAnnotatedDataDir = "syn5706203"

def updateRow(row, key="name", path="./"):
    """Remove the key from the dictionary, return it as the absolute path.
    """

    newrow = dict(row)

    filename = row[key]
    newrow.pop(key)

    return (os.path.abspath("%s/%s" % (path, filename)), newrow)

def process(row, syn, key, path, parentId):
    """Proces the row, storing and adding annotations to the file.
    
    """

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
    parser.add_argument("-p", "--parentId", type=str, help="Folder parentId to upload to.")
    parser.add_argument("datafile", type=str, help="csv file containing file name and annotations")

    args = parser.parse_args()

    syn = synapseclient.login(silent=True)

    with file(args.datafile) as f:
        reader = csv.DictReader(f)

        map(lambda x: process(x, syn, key=args.key, path=args.directory,
                              parentId=args.parentId),
            reader)

if __name__ == "__main__":
    main()
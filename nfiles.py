#!/usr/bin/env python

__author__ = 'daniela'

import glob
import argparse

def nfiles(fileroot):
    """
    Find the number of files with a specific filename.
    """

    filenames = glob.glob("*fileroot*")

    print("There are " + str(len(filenames)) + " files with root *" + str(fileroot) + "*.")
    return


def main():
    nfiles(fileroot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Find the number of files with a specific file name.")
    parser.add_argument("-f", "--froot", action="store", dest="fileroot", required=True,
                        help="Specify the fragment of the file name to search for.")

    clargs = parser.parse_args()
    fileroot = clargs.fileroot

    main()
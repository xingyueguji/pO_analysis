#!/usr/bin/env python

"""Script to replace original file names with numbered names, based on output from makeManifest pass on TeX output
    """

import re
import os
import string
import glob

def removeBraces(name):
    return name.replace('{','').replace('}','')

def checkDuplicate(name, namesList):
    """Check to see if there is a match (up to, and possibly through, the extension) for the [oldName, newName] pairs in namesList.
        Returns the match (including) or empty and sets the match[0] to name.
        """

    returnVal = ''
    for testName in namesList:
        if (name == testName[0][0:len(name)]):
            testName[0] = name
            returnVal = name
            break
    return returnVal

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--no-mv", dest="no_mv", default=False, action="store_true", help="do not execute git mv")
    parser.add_argument("texfile", type=str, help="name of .tex file")
    parser.add_argument("figfile", type=str, help="path to FigList.txt")
    args = parser.parse_args()
    with open(args.figfile,"r",encoding="utf-8") as f:
        lines = f.readlines()
        if (lines[0] == '' or lines[0]=='\n'):
            lines = lines[1:] # skip blank leading line?
        lines = [l.rstrip('\n') for l in lines]
        print(lines)
        names = [l.rsplit(' ') for l in lines] # no embedded blanks allowed!
        if not args.no_mv:
            for line in lines:
                os.system("git mv {}".format(line))
    with open(args.texfile,"r+") as f:
        # look for filenames in the TeX that did not include an extension
        # note: {filename} is used when filename contains periods
        p = re.compile("""
            \\\\includegraphics(  # pull in all included graphics files: assumed file has been trimmed of comments
            \[.*\])?            # ignore any options
            {((?:{).*(?:})|.*?)  # pull in {filename} or filename up to any extension (non-greedy): second group
            (\..*)?}     # and look for any extension .xyz: third group; full filename is second+third
            """, re.VERBOSE)
        t = f.read()
        f.seek(0)
        m = p.findall(t)
        badnames = [removeBraces(i[1])+i[2] for i in m if not(i[2])]
        # look for a match in names
        for i in badnames:
            if not(checkDuplicate(i,names)):
                print('>>> Failed to match {0}'.format(i))
        # Now we look for a name clash
        pName = re.compile('(.*/)?(.*)')
        inputNames = [pName.match(i[0]).group(2) for i in names]
        outputNames = [pName.match(i[1]).group(2) for i in names]
        ok = True
        for name in inputNames:
            if name in outputNames:
                print('>>> Name clash! Check/fix {0}'.format(name))
                ok = False
        if not(ok):
            return
        # now handle {filename} cases (when filename contains period)
        # remove extra {} first so normal renaming will work
        for filename, suffix in [(i[1],i[2]) for i in m if i[1].startswith('{') and i[1].endswith('}')]:
            filename2 = removeBraces(filename)
            t = t.replace('{'+filename+suffix+'}','{'+filename2+suffix+'}',1) # follow below pattern to avoid non-unique shorter versions
        # final renaming
        for name in names:
            if (name[0] != name[1]):
                t = t.replace('{'+name[0]+'}','{'+name[1]+'}',1) # avoid non-unique shorter versions by requiring bracketing {}
        f.write(t)
        f.truncate()

if __name__ == "__main__":
    main()

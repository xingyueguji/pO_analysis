        # To create a standalone environment
        # create a driver file for non-tdr use

import time
import logging
import shutil
import re
from pathlib import Path

class Standalone(object):

    def __init__(self, tag=None, type=None, path=None, verbosity=0):
        self._t0 = time.time()
        # set up logging
        v = {0:logging.WARNING, 1:logging.INFO, 2:logging.DEBUG}
        logging.basicConfig(level=logging.INFO, format='%(levelname)10s:\t %(funcName)12s\t %(asctime)s\t %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(v[min([2,verbosity])]) # clamp at the maximum level
        self._logger.debug('Creating Standalone instance; Starting logging')
        self._local = path
        if not tag:
            self._tag = path.stem
        else:
            self._tag = tag
        self._doctype = self.tagParse(self._tag)

    def tagParse(self, tag):  # based on noteType in cmsnotes
        # need to add CMS Note?
        p = re.compile(
            r'((?P<ntype>[ADI])N|(?P<ptype>[A-Z]{3})|(?P<btype>B2G))-(?P<yr>\d{2})-(?P<no>\d{3})')
        m = p.match(tag)
        doctype = None

        if (not m):
            raise RuntimeError('Bad tag path: {}. Expecting something like B2G-21-001 or AN-21-002'.format(tag))
        if m.group('ntype'):
            doctype = str(m.group('ntype')).lower() + 'n'
        # checking for B2G as exception to three alpha for PAS/paper
        elif m.group('ptype') or m.group('btype'):
            doctype = 'paper'
        return doctype

    def doit(self):
        driver = 'main.tex'
        with open(self._local / driver, 'w') as f:
            with open(self._local/ 'utils/general/skeleton_start.tex') as ff:
                shutil.copyfileobj(ff, f)
            with open(self._local/ 'utils/general/definitions.tex') as ff:
                shutil.copyfileobj(ff, f)
            f.writelines(['\n\\input{{{}}}\n'.format(self._tag)])
            with open(self._local/'utils/general/skeleton_end.tex') as ff:
                shutil.copyfileobj(ff, f)
        # possible document type options for the class file:
        # tdr, note, cmspaper, am, cr, in, pas, dn
        if self._doctype == 'paper':
            clsType = 'cmspaper'
        elif self._doctype == 'personal':
            clsType = 'pas' 
        else:
            clsType = self._doctype
        with open(self._local / driver, 'r+') as f:
            text = f.read()
            text = re.sub(r'\\documentclass\[11pt,twoside,a4paper,tdr\]\{cms-tdr\}', r'\\def\\svnVersion{untracked}\\def\\svnDate{\\today}\n\\documentclass[11pt,twoside,a4paper,%s]{cms-tdr}'  % clsType, text)
            f.seek(0)
            f.write(text)

        latexmk  = """
###
$pdf_mode =1;
@default_files = ('main.tex');
ensure_path( 'TEXINPUTS', './utils/general//', './figs//' );
ensure_path( 'BSTINPUTS', './utils/general//' );
###
"""
        with open(self._local / '.latexmkrc','w') as  f:
            f.write(latexmk)

        # fix older repos
            with open(self._local / (self._tag+'.tex'), 'r+') as f:
                text = f.read()
                # replace tdr-style bib with specific bib
                text = re.sub( r'\\bibliography\{auto_generated\}',r'\\bibliography{%s}' % self._tag, text)
                f.seek(0)
                f.write(text)
                f.truncate()
            


def main(argv):

    import argparse

    parser = argparse.ArgumentParser(description='TDR standalone environment creation')
    parser.add_argument('-v', '--verbosity', action='count', dest='verbose', default=False,
                        help='trace script execution: default is WARN, using -v will increase through INFO, DEBUG; going beyond that point turns on git debugging.')
    parser.add_argument('-t', '--type', action='store', dest='type', default=None,
                        help='note type. paper, pas, an, dn, etc. Usually autodetected but for pas/paper,  paper will be assumed.')
    parser.add_argument( '-x', '--tag' , action='store',  dest='tag', default=None, help='the document tag. Should be autodetected')

    opts = parser.parse_args()
    thisfile = sys.argv[0]
    local = Path(thisfile).resolve().parents[2]

    o = Standalone(tag=opts.tag, type= opts.type, path= local)
    o.doit()

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])

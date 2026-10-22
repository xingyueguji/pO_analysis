import re
from pathlib import Path
import itertools
from collections import namedtuple
import logging
import os
import subprocess
import shutil




class getFigInfo(object):
    """
    Extract information on files included via \includegraphics from the TeX log file.
    
    Creates a local named ntuple of the figures and their generated names, eg, Figure_001-a
    """

    def __init__(self, texlog: str):
        """
            :param texlog: the TeX log file names
        """

        self.texlog = texlog
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.WARNING)
        self.thumnails = True
        self.special_files = re.compile('BigDraft|cms_draft_.*|cms_pas\.*|ORCIDiD_iconvector|CMS-bw-logo|cernlogo') 

        # find the convert (ImageMagick) app
        # use default location on Windows; expect in PATH for Unix
        self.my_env = dict()
        if os.name == 'nt':
            imdir = Path(os.getenv("LOCALAPPDATA"))/"Programs"/"ImageMagick-7.1.0-Q16-HDRI"
            gsdir = Path(os.getenv("LOCALAPPDATA"))/"Programs/gs9.55.0/bin"
            if not ((imdir /"magick.exe").exists() and (gsdir/"gswin64.exe").exists()):
                print("ImageMagick or Ghostscript not found at %s and %s", imdir, gsdir)
            self.my_env["PATH"] = str(imdir) + os.pathsep + str(gsdir)
            self.convert = str(imdir/"magick.exe")
        else:
            self.convert = "convert" # assume in default path except for Windows
            if not shutil.which(convert):
                self.logger.warning("The convert program is not found in the default path")

    def fig_filter(self, s: str)->bool:
        """
        filter out known ignorable file names
        :param s: file name to check
        :returns: False for ignoreable
        """
        # Here we just look for one of the special files. 
        return not self.special_files.match(s)
    

    def subfig_grouper(self, tup: namedtuple)->str:
        """
        Returns an ordering string from a fig_info tuple

        :param tup: an instance of the namedtuple of figure information
        :returns: a string to be used in ordering /grouping such tuples
        """

        # parse the first field of the tuple, the figure 'number'
        # subfigures end with lowercase subfig part, e.g., 3a.
        # The figure number would be one less -> 2.
        p1 = re.compile('[a-z]') 
        p2 = re.compile('(?P<fignum>\d+)(?P<part>[a-z])?')
        if not p1.match(tup.num[-1:]): # last character not a letter?
            num = tup[0] 
        else:
            m = p2.match(tup.num)
            num = chr(int(m.group('fignum'))-1) # reduce the fig number part by one character

        return tup.app_label+num

    def numberwithin_check(self, b: str)->bool:
        """
        Check to make sure than any numberwithin isn't going to break everything
        
        :param b: the log file content
        :returns: OK if everything is basically monotonically increasing
        """

        # Exemplar: <789FIG ORCIDiD_iconvector.pdf A.0>
        fignums = itertools.pairwise(re.findall("""<789FIG\s*           # starting <789FIG<space>
            \S*\s*([A-Z]\.)?                                            # <any non-spaces: file name><space><optional [A-Z]. for appendices>
            (.*)                                                        # the number, with possible subpart letter
            >                                                           # closing bracket
            """, b,  re.VERBOSE|re.MULTILINE))

        for pair in fignums:
            if int(pair[1][1])< int(pair[0][1]):
                self.logger.warning("Internal figure numbers out of order: %s", pair)
                if pair[0][1][:-1] == pair[1][1]:
                    self.logger.warning("Looks like subfig. Probably will be okay.")
        return

    def create_thumbnail(self, graphic:Path, destination:Path):
        """
        Run a figure through ImageMagick aka convert to produce a thumbnail image
        
        :param graphic: the Path to the file
        :param destination: the Path to the receiving directory
        :return : status of convert subprocess
        """

        # ImageMagick arguments [https://imagemagick.org/script/command-line-options.php]
        # verbose: outputs diagnostic information
        # trim: remove area from edges the same color as the four corners. Ditch? Yes.
        # thumbnail x240: remove unneeded informtion and scale to a height of 240 pixels
        # define pdf:use-cropbox-true: use the PDF cropbox to trim
        cmd = [self.convert, str(graphic), "-verbose", "-thumbnail", "x240", "-define", "pdf:use-cropbox=true", str(destination/(graphic.stem+'-thumb.png')) ]

        s = subprocess.run(cmd, capture_output=True, encoding='utf8', env=self.my_env)
        if s.returncode != 0:
            print(s.stderr)

        return s.returncode        

    def generate_thumbnails(self, target_dir: Path, destination:Path):
        """
        generate all the thumbnails
        
        :param target_dir: Path to the base directory for the TeX build producing the log file
        :param destination: Path to the destination directory
        """

        for fi in self.named_figs:
            f = target_dir/fi.path
            if not f.exists():
                fig2 = f.parent/"fig"/f.name
                if not fig2.exists(): # check in the implicitly included directories
                    self.logger.warning(f"Fig file %s does not exist", f)
                else:
                    f = fig2

            shutil.copyfile(f, destination/fi.canonical_name)
            self.create_thumbnail(destination/fi.canonical_name, destination) 
            
   
    class figInfo(namedtuple('figInfo', ['num', 'path', 'name', 'fullname', 'app_label', 'canonical_name'])):
        """
        Named ntuple for figure info.

        Fields: ['num', 'path', 'name', 'fullname', 'app_label', 'canonical_name']
        """
        __slots__ = () # I don'think this will make much, if any, difference in efficiency, but it was in the example

        def set_canonical_name(self, figl: int, subfig_offset: int)->namedtuple:
            r"""
            Create the canonical name given the figure number for a named tuple 
            NOTE: make sure any \numberwithin is placed before the \appendix command in the TeX.
            
            :param subfig_offset: the subfigure index
            :returns: the updated tuple
            """

            fig_offset = 1; # usually the appropriate choice for current TeX distros; 0 earlier
            p = re.compile("""
                (?P<fnum>\d+)(?P<subfig>[a-z])? # break fignum into number plus possible subfig label (a, b, etc.)
                """,
                re.VERBOSE)

            m = p.match(self.num)

            if figl == 1:   # => no subfigures
                if not self.app_label:
                    tag = "Figure_{:0>3}".format(int(self.num)+fig_offset)
                else:
                    tag = "Figure_{}{:0>3}".format(self.app_label[:-1], int(m.group("fnum"))+fig_offset)
            else:
                if m:
                    if not m.group('subfig'): # regular figure
                        if self.app_label: # appendix?
                            tag = "Figure_{}{:0>3}-{}".format(self.app_label[:-1],  int(m.group('fnum'))+fig_offset, chr(ord('a')+subfig_offset)) # remove the '.' from app_label
                        else:
                            tag = "Figure_{:0>3}-{}".format(int(m.group('fnum'))+fig_offset, chr(ord('a')+subfig_offset))
                    else: # from subfigure environment
                        if m.group('subfig'):
                            tag = "Figure_{0:0>3}-{1}".format(m.group('digit'), m.group('subfig'))
                        else:
                            tag = "Figure_{0:0>3}-{1}".format(int(m.group('digit'))+fig_offset, chr(ord('a')+subfig_offset))
            fpath = Path(self.path).with_stem(tag)
            return self._replace(canonical_name = fpath.name)


    def print_files(self, target_dir: str):
        """
        print the map of file stems to canonical file names
        
        :param target_dir: base directory for the TeX build producing the log file
        """

        target_dir = Path(target_dir)
        for fi in self.named_figs:
            f = target_dir/fi.path
            if not f.exists:
                self.logger.warning(f"Fig file %s does not exist", f)
            print (f'{fi.path} {fi.canonical_name}')
            
            
    def parse_log(self) -> None:
        """ 
        Parse the TeX log file, looking for the information on figures and create a list of them.

        Adds the named_figs and other_figs items
        
        """

        with open(self.texlog, encoding='utf-8') as f:
            tlog = f.read()

        self.numberwithin_check(tlog) 

        simple = not re.search('<789FIG', tlog) #figures are already named and there are no <789FIG markers. Just extract names.
        if simple:
            p = re.compile('<.*?>\nFile:\s.*?\n<use .*?>', re.MULTILINE)
            figs = iter(p.findall(tlog))
            p = re.compile("""
                (?P<fignum>)    # no fignum  
                <(?P<figname>.*?),         
                .*?File:\s*
                (?P<full_figname>.*?) # another instance of the name as a check
                \s*Graphic.*?\n<use\s*           # this should skip over any PDF inclusion error messages
                (?P<figpath>.*?)    # look for figure path /n<use ...>
                >
                """,   
                re.VERBOSE|re.DOTALL)

        else:
            figs = iter(tlog.split("<789FIG ")) # cms-tdr.cls has redefined \includegraphics to emit <789FIG for each figure. The figure number is incremented by \caption.
            next(figs, None) # skip the stuff before the first entry
            p = re.compile("""
                (?P<figname>\S*)    # figure name, e.g., 'cernlogo'
                \s*
                (?P<applabel>([A-Z]\.)?)        # optional appendix name plus dot (eg, A.)
                (?P<fignum>.*?)>    # up to closing, including possible subfig part (eg, 3d) >
                .*?File:\s*
                (?P<full_figname>.*?) # another instance of the name as a check
                \s*Graphic.*?\n<use\s*           # this should skip over any PDF inclusion error messages
                (?P<figpath>.*?)    # look for figure path /n<use ...>
                >
                """,   
                re.VERBOSE|re.DOTALL)
        """
            <789FIG cernlogo 0>
            <cernlogo.pdf, id=67, 79.79813pt x 80.3pt>
            File: cernlogo.pdf Graphic file (type pdf)
            <use cernlogo.pdf>
            Package pdftex.def Info: cernlogo.pdf  used on input line 44.
            (pdftex.def)             Requested size: 28.45274pt x 28.63106pt.
            """

        try:
            fig_list = [self.figInfo(m.group("fignum"), m.group("figpath"), m.group("figname"), m.group("full_figname"), m.group('applabel'), "") for m in map(lambda z: p.match(z), figs)]
        except AttributeError:
            self.logger.warning("No figures found") # as there should always be a logo, at a minimum, no files implies busted TeX build
            self.named_figs = None
            return

        fig_list_filtered = filter(lambda z: self.fig_filter(z.name), fig_list)   # remove the known ignorable files
        if simple:
            self.named_figs = list(fig_list_filtered)
        else:
            grouped_figs = itertools.groupby(fig_list_filtered, self.subfig_grouper)    # group on the figure number
            grouped_figs = [ (yy[0], list(yy[1])) for yy in grouped_figs]     # convert to lists from iterators for reuse
            self.named_figs = [ z[1].set_canonical_name(len(g[1]), z[0]) for g in grouped_figs for z in enumerate(g[1])] # get list with formatted figure names, Figure_00...

        # as an extra, pull out a list of the special files present
        specials = filter(lambda z: not self.fig_filter(z.name), fig_list)
        grouped_specials = itertools.groupby(specials, lambda z: z.name)
        self.other_figs = [ next(z[1]) for z in grouped_specials] # should only be one for each special file, so use next to convert generator to value

        # get the other included files
        # only works if simple ==  False
        self.other_files = re.findall("""<567INP\s+# look for opening <567INP\s
            (             # with a (non-capturing) group composed of either
                (?:.*)\.tex       # a file with a .tex extension
                |           #or 
                (?:[^.]*)         # a file with no extension
                )>          # final closing >
            """, tlog, re.VERBOSE)

        return


if __name__ == '__main__':
    import sys
    import argparse

    if sys.hexversion <  0x3080000:
        sys.exit('This program requires python version 3.8 or above to run correctly. Try lxplus9 next time. Exiting...')

    parser = argparse.ArgumentParser(description="Get info on files, both figures and other TeX input files, from a TeX log")
    parser.add_argument( 'texlog', action='store', help='path to the TeX log file', default="testdata/template_tmp.log")

    opts = parser.parse_args(sys.argv[1:])
    fi = getFigInfo(opts.texlog)
    fi.parse_log()
    if fi.named_figs:
        for f in fi.named_figs:
            print(f)

        fi.print_files(Path("C:/Users/galve/source/repos/tdr/utils/test/testdata"))
        fi.generate_thumbnails(Path("C:/Users/galve/source/repos/tdr/utils/test/testdata"), Path("C:/Users/galve/source/repos/tdr/utils/test/testdata"))
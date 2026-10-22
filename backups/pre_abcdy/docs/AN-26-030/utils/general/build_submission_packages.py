from difflib import Match
from pydoc import describe
import re
import shutil
import socket
import sys
import os
import string
import subprocess
import tempfile
from pathlib import Path
import os
import time
import re
import json
import logging
import time
from collections import namedtuple

### Fill in the necessary command options with the correct values in the configuration file (see command line options for current name of file) or supply via the command line ###

"""
{
    "target": "EXO-21-010",
    "journal" : "JHEP",
    "ppn" : "2022-103",
    "has_AL" : false,
    "verbosity" : 1

} 
"""

#######################################################

def with_stem(path: Path, stem: str) -> Path:
    """
    Looking forward to the Python 3.9 Path.with_stem
    Assumes that there is already some suffix so that the stem is the next to last piece of the path.

    :param path: the Path
    :param stem: new stem
    :returns: the Path with the new stem but retaining all the rest, including the suffix.

    """
    parts = list(path.parts)  # final part is the name = stem + suffix
    parts[-1] = stem + path.suffix
    return Path('').joinpath(*parts)

def extract_balanced(text: str, delim: str ='{')-> list:
    """ Extract a delimited section of text: available opening delimiters are '{', '"', and  '<' 
    
        :param text: the searched text body
        :param delim: the delimiter to use. Default '{'
        :return [pout, delimited text]: [location in the search body past the matched delimiter, the delimited text (None, if not found)]
        """
    delims: dict = {"{":"}", '"':'"', "<":">"} # matching closing delims
    pin: int = text.find(delim) + 1
    nbraces: int = 1
    pout: int = pin
    while nbraces > 0:
        if pout > len(text): 
            return [0, None] # probably unmatched } inside TeX comment string
        if text[pout:pout + 2] == '\\' + delim: # look for escaped delim
            pout += 2
        else:
            if text[pout:pout + 2] == '\\' + delims[delim]:
                pout += 2
            else:
                if text[pout:pout + 1] == delims[delim]:
                    nbraces -= 1
                elif text[pout:pout + 1] == delim:
                    nbraces += 1
                pout += 1
    return [pout, text[pin:pout - 1]]

def extract_tex_field(field: str, text:str) -> str:
    """
    extract the value of a delimited TeX field
    :param field: like \\title
    :param text: the TeX source text
    :return: the value or None 
    """

    # seems a litte safer /quicker to look for the field rather than regular expressions
    pos = text.find(field)
    match = extract_balanced(text[pos:])

    return match[1]

def buildArchive(cmd: str, target: str, proc_args: dict) ->namedtuple:
    """
        tdr command line runner

        :param cmd: the basic command to tdr
        :param target: the identifier of the document to build [XXX-08-000]
        :param proc_args: standard arguments for subprocesses
        :return: named tuple (ResultFiles) of paths to the archive ['archive'], the PDF output ['pdf'], and the HTML check file ['html']

    """

    def check_file(test: str, text: str) -> Path:
        """
        Check validity of search return and pull value as a path
        
        :param test: the search string to run
        :param text: the string to search within
        :return: the resulting path or None
        """
        r = re.search(test, text, flags=re.MULTILINE)
        if r:
            f = r.group(1)
            if f:
                p = Path(f)
                if p.exists():
                    return p
        return None



    # build the submission package
    build = time.time()

    # get the logger
    logger = logging.getLogger(__name__)

    p = subprocess.run(cmd+['b', target], errors='ignore', **proc_args)
    logger.info('    Built.  Time = {:0.2f} s'.format(time.time()-build))
    if (p.stderr):
        logger.warning(p.stderr)
    logger.debug(p.stdout)
    appmiss = re.search(r' destination with the same identifier \(name\{appendix\.', p.stdout, flags=re.MULTILINE)
    if appmiss:
        logger.warning('The log file indicates that there is an unflagged appendix')
    zarc = check_file('^Contents of (.*?):$', p.stdout)
    pdf = check_file(r'^ PDF Output file:\s*(.*?)$', p.stdout)
    html = check_file(r'>>> HTML Check file:\s*(.*?)$', p.stdout)

    build_log_file = (Path('tdr_build.log')).resolve()
    with open(build_log_file, mode='w', encoding='utf-8') as build_log:
        build_log.write(20*'='+'STDERR'+20*'='+"\n")
        build_log.write(p.stderr)
        build_log.write(20*'='+'STDOUT'+20*'='+"\n")
        build_log.write(p.stdout)
        build_log.write(20*'='+'END STDOUT'+20*'='+"\n")

    results = namedtuple('ResultFiles', ['archive', 'pdf', 'html','log'])
    return results(zarc, pdf, html, build_log_file)

def set_logger(logfile: str, target: str, verbosity:int):
    """
    Set standard defaults for logging system, both file and console
    
    :param logfile: base name of logfile
    :param target: appended to logfile name
    :param verbosity: for console log stream, 0=> only log warnings and worse ... 2=> log at debug level.
    :returns: nothing
    """

    def_folder = Path.home()/'.BuildSubmissionLogs'
    v = {0:logging.WARNING, 1:logging.INFO, 2:logging.DEBUG}

    #
    # set up a console logger to replace the stock root handler; required if any of the handlers are to have a level less severe than the root
    #
    logger = logging.getLogger(__name__)
    logger.root.handlers = []
    #
    clogger = logging.StreamHandler(sys.stdout)
    clogger.setFormatter(logging.Formatter('> %(asctime)s\t %(message)s', datefmt="%Y-%m-%d %H:%M:%S"))
    clogger.setLevel(v[min([2,verbosity])]) # clamp at the maximum level
    #
    # set up a file logger for complete tracing  
    #      
    flogname = Path(logfile)
    if not Path(flogname).is_absolute():
        if not def_folder.exists():
            def_folder.mkdir()
        flogname = def_folder/flogname
    fstem = "_".join([flogname.stem, target,  time.strftime("%Y%m%d")])
    flogname = with_stem(flogname, fstem) # flogname.with_stem(fstem)
    if not (flogname.suffix):
        flogname = flogname.with_suffix(".log")
    #
    flogger = logging.FileHandler(filename=flogname) # default for logging is to append
    flogger.setFormatter(logging.Formatter('%(levelname)10s:\t %(funcName)12s\t %(asctime)s\t %(message)s', datefmt="%Y-%m-%d %H:%M:%S"))
    flogger.setLevel(level=logging.DEBUG)
    #
    logger.addHandler(flogger)
    logger.addHandler(clogger)


def get_metadata(texfile: str) -> tuple:
    """ 
    Python version of the genPreview subroutine in makeManifest.
    
    
    :param texfile: full path to the input TeX file used for the final tdr build pass
    :return (pdftitle, pdfauthor, title, abstract, m):

    """

    # Header for the HTML output

    preface = r"""
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
    <!DOCTYPE html>
    <html>
    <head>
    <title>MathJax TeX Test Page</title>
    <script type="text/x-mathjax-config">
    MathJax.Hub.Config({tex2jax: {inlineMath: [['\$','\$'], ['\\\\(','\\\\)']]}});
    </script>
    <script type="text/javascript"
    src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
    </script> </head> <body> <p>Below are the title and abstract after
    passing through the same formatter [MathJaX] as is used on CDS. You
    should be able to spot those instances where comments (%) or macros
    have sneaked into the abstract. Unrecognized macros in math mode
    will be in red, but those in text mode will not stand out.</p> <hr>
    """

    with open(Path(texfile).with_suffix('.tex')) as f:
        body = f.read()

    # extracting information from the hypertext command: hypersetup{...}
    hyper = extract_tex_field('\\hypersetup', body)
    if hyper:
        m = re.search(r'pdftitle\s*=\s*(.*)', hyper)
        if m:
            pdftitle = (extract_balanced(m.group(1))[1]).strip()
        m = re.search(r'pdfauthor\s*=\s*(.*)', hyper)
        if m:
            pdfauthor = (extract_balanced(m.group(1))[1]).strip()

    # pull plain TeX fields. Beware of those who comment out stale versions! As we look for the first, normally not a problem
    title = (extract_tex_field('\\title', body)).strip()
    abstract = (extract_tex_field('\\abstract', body)).strip()
    abstract = re.sub("\n",' ',abstract) # no newlines
    abstract = re.sub(r"\s{2,}",' ',abstract) # no multiple spaces
    # do a little cleanup of the abstract


    if abstract:
        m = re.search('(?<!\\\\)%', abstract)
        if m:
            print(">>> -------------------------------------------------------------------------- <<<\n")
            print(">>>!!!    The abstract appears to have an unescaped TeX comment character (%). <<<\n")
            print(">>>       No comments are allowed within the abstract!                         <<<\n")
            print(">>> -------------------------------------------------------------------------- <<<\n")
        with tempfile.NamedTemporaryFile(prefix='tdr-html-', suffix='.html', delete=False, mode='w+b') as htmlfile:
            htmlfile.write(preface.encode(encoding='utf-8'))
            htmlfile.write("\n".encode(encoding='utf-8'))
            htmlfile.write(title.encode(encoding='utf-8'))
            htmlfile.write("\n<br/><br/>\n".encode(encoding='utf-8'))
            htmlfile.write(abstract.encode(encoding='utf-8'))
            htmlfile.write('\n</body>'.encode(encoding='utf-8'))

    results = namedtuple('Metadata', ['pdftitle', 'pdfauthor', 'title','abstract','html'])    
    return results(pdftitle, pdfauthor, title, abstract, htmlfile.name)
            

def main(argv: list):
    """
    CMS submission journal/arXiv submission package builder
    
    Arguments are taken from the configuration file plus command line arguments. Command line overwrites file arguments.

    Produces one archive for arXiv, one for the journal, and (optionally) one supplement. For the supplement, the output PDF is to be uploaded, not the archive.

    Remember that the basic command passed to tdr will end up like:
    tdr --message='Submitted to some journal' --cernNo='year and number from cern preprint number, CERN-EP-[yyyy-nnn]' --style paper 
    PLUS
    1.  the journal processing command, (eg, --aps -)
    2.  the arXiv flag, (--arxiv)
    others are added for supplement, appendix, and error checking (--preview)

    See the logfile for complete output from the tdr processing.
    """

    import argparse

    start = time.time()

    journalname = {
        'plb':'Physics Letters B', 
        'prp':'Physics Reports', 
        'nima':'Nuc. Inst. and Meth. in Physics, A',
        'jhep':'the Journal of High Energy Physics', 
        'prl':'Physical Review Letters', 
        'epjc':'the European Physical Journal C', 
        'prd':'Physical Review D', 
        'prc':'Physical Review C', 
        'jinst':'the Journal of Instrumentation', 
        'prdrc':'Physical Review D (Rapid Communications)', 
        'npaha':'Nature Physics', 
        'natus':'Nature', 
        'rpph': 'Reports on Progress in Physics',
        'csbs': 'Computing and Software for Big Science'}
    journalcmd = {
        'plb': '--plb=3p,twocolumn,times', 
        'prp': '--prp=3p,twocolumn,times', 
        'nima': '--nima=3p,twocolumn,times',         
        'prl': '--aps=reprint,prl,longbibliography', 
        'jhep': '--jhep=-', 
        'epjc': '--epjc=twocolumn', 
        'prd': '--aps=reprint,prd,longbibliography', 
        'prdrc': '--aps=reprint,prd,longbibliography', 
        'prc': '--aps=reprint,prc,longbibliography', 
        'jinst': '--jhep=-', 
        'npaha': '--npaha=-', 
        'natus': '--npaha=-', 
        'rpph' : '--rpp=-',
        'csbs' : '--csbs=-'
        }

    # default option values
    def_gitUrl = 'ssh://git@gitlab.cern.ch:7999/'
    def_conffile = 'BuildSubmission.json'
    def_has_AL = False
    def_logfile = 'SubmissionBuilder'
    def_remote = True
    def_journal = 'JHEP'

    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument( 'target', action='store', default=False, help='paper identifier [XXX-08-000]')
    parser.add_argument('-v', '--verbosity', action='count', dest='verbosity', default=1, help='trace script execution: default is INFO, using additional -v will increase the level')
    parser.add_argument( '-c', '--conf', action='store', dest='conffile', default=def_conffile, help=f'configuration filename [default: {def_conffile}]')
    parser.add_argument( '-r', '--remote', action=argparse.BooleanOptionalAction, dest='use_remote', default=def_remote, help=f'Merge in json from https://gitlab.cern.ch/tdr/papers/target; default: {def_remote}. Negate with --no-remote.')
    ca = parser.add_argument_group('local configuration parameters')
    ca.add_argument('-o', '--arch_dir',  action='store', dest='arch_dir', default=None, help='Output directory for final zip/tar files and PDF. Default takes TDR default')
    ca.add_argument('--own_tdr', action='store', dest='own_tdr', default=None, nargs='?', help='use a local tdr: requires full path, /.../tdr, [default: None => directory tree of this script]')
    ca.add_argument('--own_repo', action='store', dest='own_repo', default=None, nargs='?', help="Path to local tdr/papers/target to use in place of pulling from GitLab; fault=None=> cwd")
    ca.add_argument('--logfile', action='store', dest='logfile', default=def_logfile, help=f'file name stem for logger output; will have target and timestamp added [default: {def_logfile}]')
    ca.add_argument('--gitUrl', action='store', dest='git_url', default=def_gitUrl, help=f'URL to use for git outside of CERN. At CERN use Kerberos authentication. [default outside CERN: {def_gitUrl}]' )
    pa = parser.add_argument_group('paper parameters (ppn is required)')
    pa.add_argument('--journal', action='store', dest='journal', default=def_journal, help=f'target journal [jhep, prl, prd, prc, epjc, plb, prp, etc.]; default: {def_journal}' )
    pa.add_argument('--has_app', action=argparse.BooleanOptionalAction, dest='has_app', default=False, help='Has an appendix? [default: False] ' )
    pa.add_argument('--ppn', action='store', dest='ppn',  help='CERN preprint number w/o leading CERN-EP-')
    pa.add_argument('--has_supp', action=argparse.BooleanOptionalAction, dest='has_supp', default=False, help='Has a supplement? [default: False]' )
    pa.add_argument('--doi', action='store', dest ='doi', default=False, help='DOI for published papers w/o leading dx.doi.org, eg [10.1007/JHEP05(2022)014]')
    pa.add_argument('--has_AL', action=argparse.BooleanOptionalAction, dest='has_AL', default=def_has_AL, help=f'Should paper include author list/acks? [default: {def_has_AL}. Negate with --no-has_AL. Hopefully a temporary measure]')
    pa.add_argument('--suppJFormat', action=argparse.BooleanOptionalAction, dest='suppJFormat', default=True, help='Build supplement in journal format [False=> CMS format; default: True. Negate with --no-suppJFormat.]')

    opts = parser.parse_args(argv) # parse them all to get a possible configuration file

    # combining a config file and command line arguments. Ref: https://medium.com/swlh/efficient-python-user-interfaces-combining-json-with-argparse-8bff716f31e4
    if Path(opts.conffile).exists(): # and not opts.use_remote:
        opts.conffile = Path(opts.conffile).absolute()
        with open(opts.conffile) as f:
            config = json.load(f)
        vars(opts).update(config) # the vars operator will unpack the opts Namespace, allowing us to update the values from those in config
        parser.parse_args() # override from command line arguments by reparsing

    # set up logging
    base_log_level = logging.DEBUG # set to fairly permissive as base
    logging.basicConfig(level=base_log_level)
    logger = logging.getLogger(__name__)    
    if opts.logfile:
        set_logger(opts.logfile, opts.target, opts.verbosity)
    else:
        v = {0:logging.WARNING, 1:logging.INFO, 2:logging.DEBUG}
        logger.setLevel(v[min([2,opts.verbosity])])
    logger.info(f'Starting build process for {opts.target}')

   
    # evaluate any relative paths in the options
    if not opts.use_remote:
        if opts.own_repo:
            opts.own_repo = Path(opts.own_repo).resolve()
        elif opts.own_repo is None:
            opts.own_repo = Path.cwd().resolve()
        if opts.own_tdr:
            opts.own_tdr = Path(opts.own_tdr).resolve()
        elif opts.own_tdr is None:
            opts.own_tdr = (Path(__file__).parents[1]/'tdr').resolve()

    # get a clean copy of the paper into a working temp area
    wd = Path.cwd()

    twd = Path(tempfile.mkdtemp(prefix='gitTemp-')).absolute()
    os.chdir(twd)
    logger.info(">>> Building submission packages for {} in temporary directory {}...".format(opts.target, twd))
    try: # try, but not too hard
        logger.info(">>>          with logfile {}".format(logger.handlers[0].baseFilename))
    except:
        pass
    host = socket.getfqdn()
    my_env = os.environ.copy()
    proc_args =  {'stdout':subprocess.PIPE, 'stderr':subprocess.PIPE, 'check':True, 'encoding':'utf-8', 'env':my_env} # capture STDOUT and STDERR
    fetch = time.time()

    # get the repo for the target paper and for utils
    url = 'https://:@gitlab.cern.ch:8443/'
    if host.endswith('.cern.ch') and host.startswith('lxplus7'):    
        my_env["PATH"] = "/opt/rh/rh-git29/root/usr/bin:" + my_env["PATH"]
        my_env["PERL5LIB"] = "/opt/rh/rh-git29/root/usr/share/perl5/vendor_perl"
        my_env["LD_LIBRARY_PATH"] = "/opt/rh/httpd24/root/usr/lib64"
    elif host.endswith('.cern.ch') and host.startswith('lxplus'): 
        pass
    else: # assume everything (TeX, perl, etc) is correctly installed and we are using ssh in place of https
        url = opts.git_url
        proc_args.pop('env')
    if url[-1:] != '/': # make trailing / optional by adding it here
        url += '/'
    logger.info('Using URL %s', url)

    try:
        # own copy of the target repo?
        if not opts.own_repo:
            p = subprocess.run(["git", "clone", url+"tdr/papers/"+opts.target], **proc_args)
            if (p.stderr and (p.stderr != f"Cloning into '{opts.target}'...\n" and p.stderr != f"Cloning into '{opts.target}'...\nwarning: redirecting to https://gitlab.cern.ch:8443/tdr/papers/{opts.target}.git/\n")):
                logger.warning(p.stderr)
            logger.debug(p.stdout)
        else:
            repo = Path(opts.own_repo) #/ opts.target
            logger.info("Using local copy of target repo, %s", repo.resolve())
            if repo.exists() and repo.is_dir():
                shutil.copytree(repo, twd / opts.target)

        # own copy of tdr/utils?
        if not opts.own_tdr:
            p = subprocess.run(['git', '-C', str((Path(__file__).resolve()).parent), 'rev-parse', '--abbrev-ref', 'HEAD'], capture_output=True)
            if p.returncode == 0:
                branch = p.stdout.decode(encoding='utf-8').strip()
                p = subprocess.run(["git", "clone", "-b", branch, url+"tdr/utils"], **proc_args)  
            else:
                p = subprocess.run(["git", "clone", "-b", 'feature/ziploc', url+"tdr/utils"], **proc_args)  # this now picks up the development version
            if (p.stderr and (p.stderr != "Cloning into 'utils'...\n" and p.stderr != "Cloning into 'utils'...\nwarning: redirecting to https://gitlab.cern.ch:8443/tdr/utils.git/\n")):   
                logger.warning(p.stderr) 
            logger.debug(p.stdout)
            tdr = twd / 'utils/tdr' # default location in repo space
        else:
            tdr = Path(opts.own_tdr)
        if not (tdr.is_file and os.access(tdr, os.X_OK)):  # check that it actually exists and is executable
                raise NameError(' TDR executable {s} no good'.format(tdr))
        logger.info('    Target and tdr repos ready. Time = {:0.2f} s'.format(time.time()-fetch))
    except subprocess.CalledProcessError as e:
        logger.fatal(f"Problem with Git. Error code {e.returncode}:\n\n{e.stderr}")
        sys.exit() # just exit rather than re-throw

    if opts.use_remote:
        # fold in any config commands from the repo itself.
        if (Path(opts.target)/def_conffile).exists():
            opts.conffile = Path(opts.conffile).absolute()
            with open(Path(opts.target)/def_conffile) as f:
                config = json.load(f)
            vars(opts).update(config)
            parser.parse_args()
        else:
            logger.warning('No configuration file found in remote repo')


    #
    # We can now check all the remaining options
    #

    # this happens enough that we explicitly fix it
    opts.journal = (opts.journal).lower()

    logger.debug("Arguments:\n"+"\n".join("{0}:\t{1}".format(k,v) for k,v in vars(opts).items()))

    if opts.doi:
        message = 'Published in '+journalname[opts.journal]+' as \\href{http://dx.doi.org/'+opts.doi+'}{\\doi{'+opts.doi+'}.}'
    else:
        message = 'Submitted to '+journalname[opts.journal]

    if opts.ppn:
        m = re.match(r'\d{4}-\d{3}', opts.ppn)
        if not m:
            logger.fatal('Incorrect CERN preprint number, %s', opts.ppn)
    else:
        logger.fatal('Missing CERN prepring number')

    # standard build arguments
    cmdbase = ['perl', tdr,  '--style', 'paper', '--message', message, '--cernNo', opts.ppn]
    if not opts.has_AL:
        cmdbase += ['--nouseAL']
    if opts.arch_dir:
        cmdbase += ['--arc_location', opts.arch_dir]


    # Build the archives
    os.chdir(twd / opts.target)

    # prequel: get the metadata from the TeX target
    meta = get_metadata(opts.target)
    logger.info(f">>> HTML title+abstract check file:\t{meta.html}")


    # first pass: journal

    logger.info('>>> Building journal submission package')
    cmd = cmdbase
    bd = twd / 'journal'
    os.mkdir(bd)
    if opts.has_app:
        cmd = cmdbase + ['--appendix']   # this is just for JHEP submissions   
    rj = buildArchive(cmd + ['--arxiv', journalcmd[opts.journal], '--temp_dir', str(bd)], opts.target, proc_args)


    # for arXiv
    if opts.has_app or opts.has_supp: # for the arXiv version a supplement means that there is a corresponding appendix
        cmd = cmdbase + ['--appendix']
    else:
        cmd = cmdbase
    bd = twd / 'arxiv'
    os.mkdir(bd)
    logger.info('>>> Building arXiv submission package')
    ra = buildArchive(cmd +  ['--arxiv', '--temp_dir', str(bd)], opts.target, proc_args) 

    # if separate supplement... supplement must be a separate tdr-style document with the name target_supp
    if opts.has_supp:
        if opts.suppJFormat:
            cmd = cmdbase + [journalcmd[opts.journal]]
        logger.info('>>> Building the supplement using {}'.format(Path.cwd()))
        rs = buildArchive(cmd + ['--preflight', '--supplement', '--no-draft'], opts.target+'_supp', proc_args) 

    logger.info('>>> Done with all. Total elapsed time = {:0.2f} s'.format(time.time()-start))
    logger.info('>>> Journal submission: {}'.format(rj.archive))
    logger.info('>>> ArXiv submission: {}'.format(ra.archive))
    if opts.has_supp:
        logger.info('>>> Supplement PDF: {}'.format(rs.pdf))
    logger.debug(f"\n\nSubmitted to {journalname[opts.journal]}. All figures and tables can be found at http://cms-results.web.cern.ch/cms-results/public-results/publications/{opts.target} (CMS Public Pages).\n\nCMS-{opts.target}, CERN-EP-{opts.ppn}\n")

    os.chdir(wd)
    print(12*'-',"\n>>>Snippets for arXiv metadata items\n")
    for k,v in meta._asdict().items():
        print(f'{k:>12}:\t{v}\n')
    print(f"\n>>>\nSubmitted to {journalname[opts.journal]}. All figures and tables can be found at http://cms-results.web.cern.ch/cms-results/public-results/publications/{opts.target} (CMS Public Pages).\n\nCMS-{opts.target}, CERN-EP-{opts.ppn}\n")
    print(12*'-')
    logging.shutdown()

        


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])

import requests
from pathlib import Path
import shutil
from io import BytesIO
import sys
import os



def getCollab(argv):

    """
    Get the collaboration author list from the icms server (see the url below)

    :arg filename: name of the base TeX file of the paper, w or w/o _temp
    :arg type: tex|xml, type of author list to fetch
    :arg tex_dir: the working TeX directory
    :arg dest_dir: destination directory. None implies no write, i.e., fetch as a test.
    :returns: 0 [error] or the name of the fetched file
    """

    # this is from the transition from perl
    (filename, tagtype, tex_dir, dest_dir) = argv[0:4]

    # this is for `temporary` debugging
    verbose = True

    tagfile = Path(filename).stem #remove any parent directories or suffix (.tex or .xml) for fetch
    if (tagfile.endswith('_temp')):
        tagfile = tagfile[:-5]



    if (tagtype.upper() == 'XML'):
        tagfile +=  '-authorlistN.xml' # now only use the N format XML file
    else:
        tagfile += '-authorlist.tex'


    # Check for existing version in TeX working directory

    dir = Path(tex_dir)
    if (dir/tagfile).exists():
        shutil.copyfile(dir/tagfile, Path(dest_dir)/tagfile)
        print(">>>  Used local copy of author list file: ", tagfile)
        return str(Path(dest_dir)/tagfile)

    # Fetch from server

    url = 'https://icms.cern.ch/tools-api/api/alFile'
    form = {'fName': tagfile}

    certs = Path(__file__).parent/"CERN-bundle.pem" # certificate bundle should be co-located with this file
    if not certs.exists:
        print(f'>>> Error locating CERN certificate bundle, {certs}')
        return False
    else:
        r = requests.post(url, data=form, verify=certs) # icms.cern.ch will fail verification w/o explicit cert outside CERN
        if not r.status_code==requests.codes['ok']:
            print(f'>>> Error fetching authorlist {tagfile}')
            print(r.reason)
            return 0 # for perl
        elif not dest_dir is None:
            with open(Path(dest_dir)/tagfile, 'wb') as f:
                f.write((BytesIO(r.content)).getbuffer())
            if verbose:
                print(20*'>',f"\n>>> fetched {tagfile} from the server <<<\n",20*'<',sep='')
            return str(Path(dest_dir)/tagfile)
        else:
            return True # ok fetch, but no file written


if __name__ == '__main__':
    import sys
    sys.exit(not getCollab(sys.argv[1:]))

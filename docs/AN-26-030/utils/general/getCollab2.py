
import requests
from pathlib import Path
import shutil
from io import BytesIO
import sys
import os
from hashlib import sha256


def getCollab2(argv):

    """
    Get the collaboration author list from the icms server (see the url(s) below)

    :arg filename: name of the base TeX file of the paper
    :arg tex_dir: the working TeX directory
    :arg dest_dir: destination directory. None implies no write, i.e., fetch as a test.
    :returns: 0 [error] or the name of the fetched file
    """

    # this is from the transition from perl
    (filename, tex_dir, dest_dir) = argv[0:3]
    

    # this is for `temporary` debugging
    verbose = True


    # Check for existing version in TeX working directory
    dir = Path(tex_dir)
    if (dir/filename).exists():
        shutil.copyfile(dir/filename, Path(dest_dir)/filename)
        print(">>>  Used local copy of author list file: ", filename)
        return str(Path(dest_dir)/filename)
    else:
        # Fetch from server [thanks to Andreas P.]
        token = sha256( str.encode( f'-icms-al-file-{filename}' ) ).hexdigest()

        url = f'https://pfeiffer.web.cern.ch/cgi-bin/getFile.py?fileName={filename}'
        params = { "token": token }
        headers = { "Accept": "application/text" }
        request = requests.get(url)
        
        r = requests.post(url=url, data=params, headers=headers)
        if not r.status_code==requests.codes['ok']:
            print(f'>>> Error fetching authorlist {filename}')
            print(r.reason)
            return False
        elif not dest_dir is None:
            with open(Path(dest_dir)/filename, 'wb') as f:
                f.write((BytesIO(r.content)).getbuffer())
                if verbose:
                    print(20*'>',f"\n>>> fetched {filename} from the server <<<\n",20*'<',sep='')
                return str(Path(dest_dir)/filename)
        




if __name__ == '__main__':
    import sys
    sys.exit(not getCollab2(sys.argv[1:]))

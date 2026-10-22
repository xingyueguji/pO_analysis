from lxml import etree
import logging

class Manifest(object):
    """
    Builds a manifest for a CDS upload of a CMS document
    """
    def __init__(self):
        self.logger = logging.getLogger(__name__)

        self.tag = 'XXX-08-000'

    def buildTree(self):
        """
        Build the XML tree describing the document
        """

        root = etree.Element('record')
        root.append(etree.Element('037'))
        root.append(tag="001", ind1=" ", ind2=" ")
        s = etree.tostring(root)
        pass

if __name__ == '__main__':
    import sys
    import argparse

    if sys.hexversion <  0x3080000:
        sys.exit('This program requires python version 3.8 or above to run correctly. Try lxplus9 next time. Exiting...')

    parser = argparse.ArgumentParser(description="Create a CDS manifest in XML format")
  
    manifest = Manifest()
    manifest.buildTree()
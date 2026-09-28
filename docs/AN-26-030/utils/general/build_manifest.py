from lxml import etree
import logging
import tempfile
from pathlib import Path


class Manifest(object):
    """
    Builds a manifest for a CDS upload of a CMS document
    """
    def __init__(self):
        self.logger = logging.getLogger(__name__)

        #FIXME: supply and check these
        #$baseDir, $outDir, $title, $updateRecord, $xml_style, $gitYear, $dataNotMC, $abstract, $contactAddress, $thumbnails, $artTyp
        self.tag = 'XXX-08-000'
        self.doc_type = 'paper'
        self.zip_dir = Path(tempfile.mkdtemp(prefix='xxxTemp-')).absolute()
        self.update_record = "" # record number to update, example: 2839924
        self.gitYear = "2022"
        self.dataNotMC = True
        self.thumbnails = True
        self.title = 'Test document title'
        self.abstract = 'test document abstract, La la la la...'
        self.contact_address = "cms-publication-committee-chair@cern.ch"
        self.art_type = 0
        self.root = etree.Element('record')



    def buildTree(self):
        """
        Build the XML tree describing the document

        See https://cds.cern.ch/help/admin/howto-marc for information on the format

        """

        def insert_datafield(value:str, tag: str, code:str='a', ind1:str= " ", ind2:str=" ", record:etree.Element=None)->etree.Element:
            """
            Create a CDS datafield record with a subfield child and attach to the passed element

            A typical datafield set will look like
            <datafield ind1=" " ind2=" " tag="110"><subfield code="a">CMS Collaboration</subfield></datafield>
            
            :param value: the 'text' value for the subfield child
            :param tag: the 'tag' for the datafield element
            :param code: the 'code' attribute for the child element. A value of None suppresses the child element.
            :param ind1: value for ind1. Default " "
            :param ind2: value for ind2. Default " "
            :param record: record to append to. Defaults to self.root
            returns: the datafield record
            """

            if not record:
                record = self.root

            r = etree.Element('datafield', ind1=ind1, ind2=ind2, tag=tag)
            if code:
                e = etree.Element('subfield', code=code)
                e.text = value
                r.append(e)

            record.append(r)
            return r


      
        #Note is a special article type: goes into a 980_a
        art_types =  ("Particle Physics - Experiment", "Nuclear Physics - Experiment", "Detectors and Experimental Techniques", "Note")

        # insert the original CDS number if this is to be an update
        if self.update_record: 
            e = etree.Element("controlfield",  tag="001", ind1=" ", ind2=" ")
            e.text = self.update_record
            self.root.append(e)
        # CMS note number (MARC Source of Acquisition)
        insert_datafield("CMS-{}-{}".format((self.doc_type).upper(), self.tag), '037') 
        # Title (MARC Title Statement)
        insert_datafield(self.title, '245')
        # Other standard MARC data as will be implemented by CDS (as per message 2013-03-05)
        
        # collaboration, both as collaboration and as author
        insert_datafield("CMS Collaboration", '710', 'g')
        insert_datafield("CMS Collaboration", '110', 'a')

        # Accelerator field (as per Annette Holtkamp request of 2013-03-28
        r = insert_datafield("CERN LHC", tag='693')
        e = etree.Element('subfield', code='e')
        e.text = "CMS"
        r.append(e)

        # Year (MARC Publication)
        insert_datafield(self.gitYear, tag='260', code='c')

        # MARC/CERN pubdata
        r = insert_datafield("Geneva", tag="269")
        e = etree.Element('subfield', code='b')
        e.text = "CERN"
        r.append(e)
        r = insert_datafield("CERN", tag="690", ind1="C")

        # Data/MC indicator: tag 653, ind1:1, subfield 9: CMS, subfield a: Data/Monte-Carlo: from J-Y LM, 2010/02/10
        r = insert_datafield("CMS", tag="653", code="9", ind1="1")
        if self.dataNotMC:
            dMC = 'Data'
        else:
            dMC = 'Monte-Carlo'
        e = etree.Element('subfield', code='a')
        e.text = dMC
        r.append(e)

        # Abstract
        insert_datafield(self.abstract, tag="520")

        # MARC Electronic Location and Access: ind1 = 4 =>HTTP, subfield y=> Link text
        r = insert_datafield("http://cms.cern.ch/iCMS/analysisadmin/cadi?ancode={}".format(self.tag), tag="856", ind1="4", code="u")
        e = etree.Element('subfield', code='y')
        e.text = "Additional information for the analysis (restricted access)"
        r.append(e)

        # 980 tag: (MARC Equivalence or Cross-Reference Series Personal Name/Title)"
        if self.doc_type == 'paper':
            insert_datafield("CMS_PAPERS", tag="980")
        else:
            insert_datafield("NOTE", tag="980")
            insert_datafield("CMS-PHYSICS-ANALYSIS-SUMMARIES", tag="980")

        # 859:8560 tags: (email:  subfield f)
        r = insert_datafield(self.contact_address, tag="859", code="f")
  
        # 65017 tag: CERN and document type. Options are "Particle Physics - Experiment", "Nuclear Physics - Experiment", and "Detectors and Experimental Techniques"
        if self.art_type == 2: # detectors, etc.
            insert_datafield(art_types[self.art_type], tag= "980")
        else:
            r = insert_datafield("SzGeCERN", tag="650", code="2", ind1="1", ind2="7")
            e = etree.Element('subfield', code="a")
            e.text = art_types[self.art_type]
            r.append(e)

        s = etree.tostring(self.root, pretty_print=True)
        pass

if __name__ == '__main__':
    import sys
    import argparse

    if sys.hexversion <  0x3080000:
        sys.exit('This program requires python version 3.8 or above to run correctly. Try lxplus9 next time. Exiting...')

    parser = argparse.ArgumentParser(description="Create a CDS manifest in XML format")
  
    manifest = Manifest()
    manifest.buildTree()
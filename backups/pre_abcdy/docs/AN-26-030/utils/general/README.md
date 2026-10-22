# Guide to the general folder

This folder contains items of general utility in creating and manipulating CMS documents.

## files for creating documents

|||
|---|---|
|cms-tdr.cls | LaTeX class file for all CMS TeX documents|
|cms_unsrt.bst | BibTeX style file for TeX documents (same as lucas_unsrt.bst, which is the old name)|
|ptdr_definitions.sty | CMS TeX macros|
|BigDraft.pdf | background watermark image for draft version|
|cernlogo.pdf | CERN logo|
|CMS-bw-logo.pdf | CMS logo (in black and white)|
|skeleton_*.tex | required beginning and ending files included with document body|
|cms_*.pdf | header for document type *|
|AT/LHCb/TOTEM | files for joint submissions|

## LaTeX style files

These are required style files not usually found in TeX distros

|||
|---|---|
|heppennames2.sty | defines PENNAMES2 macros for elementary particle names|
|hepparticles.sty | used by heppenames2|
|hypernat.sty | hyperlinks for natbib|
|pdfdraftcopy.sty | creates watermark for PDF draft versions|
|subdepth.sty | fine tweaks for hepparticles.sty super/subscripts|
|topcapt.sty | for captions above tables|
|ulem.sty | enhanced underlining|

## journal style files

We use the standard CMS style for JHEP. The APS style files are readily available in the RevTeX package.

|||
|---|---|
|iop*.clo/cls | IOP style files|
|sv*.clo/sv*.clo/cls | Springer (EPJC) style files|
|spphys.bst | Springer BibTeX style file|

## Bib files

|||
|---|---|
|CMSPapers.bib | all the published CMS papers in our style|
|gen.bib | a list of often-used references|
|higgs.bib | Higgs history papers|
|pasBib-tech.bib | PASs as tech reports|

## example files/documentation

|||
|---|---|
|myMacro.C | for generating figures from ROOT|
|histo.root | figures from ROOT|
|CMS_lumi.* | lumi for ROOT|
|tdrstyle.C | more|
|notes_for_authors.* | extended documentation on TDR|

## TDR internal use files

|||
|---|---|
|getCollab.pl |  fetches the official authorlist|
|makeManifest.pl | (1) generates a manifest for CDS uploads, (2) runs checks|
|cleanRefs.py | checks the bib file for errors and rewrites for journal submissions|

## utilities
|||
|---|---|
|find-8bit.py | checks for Unicode characters outside the standard ASCII character set|
|matchTeXCommands.py | scans a definition file for the list of defintions actually used|
|renameFigures.py | replace original file names with numbered names, based on output from makeManifest|
|tdrDiff.py | generates a latexdiff PDF of two different revisions of a tdr-style document from Git|
|sortbib.pl | sort bib file on key|
|standalone.py | create a standalone environment for a document|

## internal use

|||
|---|---|
|get-ptdr-defs.pl | formats the standard definitions for inclusion in documentation|

## legacy files

|||
|---|---|
|pennames.sty | the original|
|pennames-*.sty | extensions for the original PENNAMES|
|placeins.sty | used only for joint ATLAS submissions (TOP-17-006 )|

## obsolete: to be removed

|||
|---|---|
|patch_html2.pl | from the original physics TDRs|
|figstrip3.pl | ditto|
|pas-bib.py | generate bib file of PAS or papers from CDS. superseded by doclist|
|pasBib.bib | PASs (not recently updated)|

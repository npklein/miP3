                        miP3.py Version 2

Example to run
--------------
After installing the dependencies (see below) an example can be run with

	python2.7 -p exampleData/TAIR10_pep_20101214_subset -i exampleData/example_tfs.fasta -f exampleData/Pfam.txt -o exampleData/test.tsv -b /path/to/blast/dir/bin/ -s 200 -a 550 -e 0.01 -z 0.5 -x 1

Installation
-----------
Unpack the zipped file and install the requirements listed in "Dependencies". After installations of the requirements it can be run from the command line with Python. 

What is it?
-----------
miP3.py is a python program that predicts microProteins from a sequenced genome. It can also be used to find similar proteins lacking any domain specified by the user. It has to be run from the command line.

Authors
-------
Niek de Klein, Sue Rhee, Enrico Magnani

Contributors
------------
Michael Banf

Contact
------
niekdeklein@gmail.com

Dependencies
------------
Internet connection
Python 2.7.x - obtainable from http://www.python.org/download/releases/2.7/ 
Biopython - obtainable from http://biopython.org/wiki/Download 
BLAST+ - obtainable from
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/. Version 2.2.29 is the version the program was tested with. For the latest version of BLAST go to ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. NOTE: The BLAST+ program does not work when there are spaces in the filepath (e.g. if it is installed in C:\Program Files\ it will not work).
SOAPpy - Read download instructions from http://www.aheil.de/2013/08/17/soap-with-python/


System requirements and runtime
-------------------------------
miP3.py runs on all operating systems with the Dependencies installed.
Running miP3.py with all transcription factors of Arabidopsis thaliana (2296 proteins) against the complete proteome of Arabidopsis thaliana found in TAIR (35386 proteins) using BLAST 2.2.29+ on a 64-bit Windows 7 Enerprise with Intel Core i5-3570 CPU @ 3.40GHz with 8 GB of RAM takes 5 hours and 18 minutes. 


How to run
----------
miP3 can be run from the command line using Python. 

    python2.7 miP3.py

Or on Windows with default installation location

    C:\Python\python2.7.exe miP3.py

There are a number of required arguments. For more information on each argument, run

    pyhon2.7 miP3.py -h

The mandatory arguments are:

    -p: A fasta file containing all proteins 
    -i: A fasta file containing proteins of interest 
    -f: Text file with unwanted Pfam domains 
    -o: The output file name 
    -b  Folder that contains makeblastd, blastp and rpsblast

And the optional arguments are:

    -s: Max size of small proteins. It is more difficult to find homologs of small proteins, so proteins below this size have a slightly different method for homolog detection.  
    -a: Max size of all proteins to search with. Sometimes you are only interested in proteins of up to a certain size. 
    -e: E-value to use when searching for homologues with all proteins. Can be set low. 
    -z: E-value to use when searching for homologues with small proteins. 
        Because it is more difficult to find homologues with small proteins, 
        use a higher e-value. There is an extra check involved to prevent 
        non-homologs from being included. 
    -x E-value to use when checking if reblasting the results from
       finding homologs of small proteins gives transcription factors as
       top results. 


All necessary files to run miP3.py to find microProteins in Arabidopsis are located in the miP3_version_2 folder except for the ncbi_blast_2.2.29+ BLAST folder. To search for miPs in Arabidopsis do: 

    python2.7 miP3.py -p TAIR10_pep_20101214 -i arabidopsis_transcription_factors.fasta -f Pfam.txt -o miP_output.csv -b ncbi_blast_2.2.29+/bin/
    
The result is written as a tab delimited file. The first column is the name of the predicted miP. The second column contains the transcription factors that are homologs of the predicted miP. The third column contain s the domains that the predicted miP contains. The final column is the length of the predicted miP.


New in this version
----------------
Version 1 of miP3 was developed for identification of microProteins in Arabidopsis thaliana and other organisms [1]. In version 1, a lot of dependencies had to be installed locally to be able to run the program. For version 2 only Python, biopython and BLAST+ have to be installed. Additionally, the code has been cleaned up and is easier to maintain and small improvements in the implementation have been made.

Changes:
- Uses a newer BLAST+ version (2.2.29)
- Uses a newer InterproScan version (5)
- Uses the web service of InterproScan instead of the local version
- Does not do a separate local Pfam search
- Putative miPs are no longer filtered out because they are in the proteins of interest file, in case a miP is misclassified.
- Uses different default thresholds based on more thorough performance testing
- BLAST against all proteins < 550a.a. instead of against all proteins
- Reblast against a database made from proteins of interest instead of all proteins
- Domains to use as filter can now be IPR IDs instead of domain names 
- Code is reorganized and cleaned up 

1. Magnani, Enrico, Niek de Klein, Hye-In Nam, Jung-Gun Kim, Kimberly Pham, Elisa Fiume, Mary Beth Mudgett, and Seung Yon Rhee. "A comprehensive analysis of microProteins reveals their potentially widespread mechanism of transcriptional regulation." Plant Physiology (2014): pp-114.


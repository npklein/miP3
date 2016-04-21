# Copyright (c) 2014 - Niek de Klein, Enrico Magnani, Seung Y. Rhee
#
#     This file is part of microProtein Prediction Program (miP3) version 2.
#
#     microProtein Prediction Program (miP3) version 2 is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     microProtein Prediction Program (miP3) version 2 is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with miP3.py.  If not, see <http://www.gnu.org/licenses/>.")

import shutil
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbirpsblastCommandline
from Bio.Blast import NCBIXML
import re
import os
import subprocess
import ntpath


def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
    
def makeBLASTdb(fasta_file, database_name, blast_path): 
    """
    Make a BLAST database from <proteins_file> located at <blast_path> with the name <database_name>
    
    Args:
        proteins_file (str):	 Path to a fasta file that has to be used to make the database
        database_name (str):	 Name for the database
        blast_path (str):		 path where the makeblastdb program is located
        
    Returns:
        Nothing
    """
    
    try:
        # try to run it with location of blast program
        subprocess.call([blast_path+'/makeblastdb', '-in',  fasta_file, '-dbtype', 'prot', '-out', database_name])
        return
    except OSError as e:
        print(e)

    # else, try to copy it to current folder and run it directly
    shutil.copy2(blast_path+'makeblastdb','makeblastdb')
    subprocess.call(['./makeblastdb', '-in',  fasta_file, '-dbtype', 'prot', '-out', database_name])
    return 

def blastp(blast_file, blast_db, evalue, blast_path):
    """
    Run a BLASTP search with blast_file against blast_db
    
    Args:
        blast_file (str):	  Path to fasta file used to BLAST against blast_db
        blast_db (str):		  Name of the database to BLAST against 
        evalue (float):		  Evalue to use as threshold
        blast_path (str):	  Path to the blastp program
    
    Returns:
        An iterable of blast records as returned by NCBIXML.parse
    """
    def cline():
        blastp_cline()
        result_handle = open("blastpOutput.xml")
        blast_records = NCBIXML.parse(result_handle)
        return blast_records
    print('starting blast, may take a while')
    try:
        # try to run it with location of blast program
        blastp_cline = NcbiblastpCommandline(blast_path+'/blastp', query=blast_file, db=blast_db, evalue=evalue,
                                                                      outfmt=5, out="blastpOutput.xml")
        return cline()
    except OSError as e:
        print(e)

    # else, try to copy it to current folder and run it directly
    shutil.copy2(blast_path+'blastp','blastp')
    blastp_cline = NcbiblastpCommandline('./blastp', query=blast_file, db=blast_db, evalue=evalue,
                                                                      outfmt=5, out="blastpOutput.xml")
    return cline()

def rpsblast(blast_file, rpsblast_db, blast_path, evalue):
    """
    Run a rpsBLAST search with blast_file rpsblast_db
    
    Args:
        blast_file (str):	  Path to fasta file used to BLAST against rpsblast_db
        rpsblast_db (str):		 Name of the database to BLAST against 
        blast_path (str):	  Path to the blastp program
        evalue (float):		  Evalue to use as threshold 
    
    Returns:
        An iterable of blast records as returned by NCBIXML.parse
    """
    def cline():
        rpsblast_cline()
        result_handle = open("rpsblastOutput.xml")
        blast_records = NCBIXML.parse(result_handle)
        return blast_records
    
    try:
        # first try to run it with location of blast program
        rpsblast_cline = NcbirpsblastCommandline(blast_path+'/rpsblast', query=blast_file, db=rpsblast_db, evalue=0.00000001,outfmt=5, out="rpsblastOutput.xml")
        return rpsblast_cline()
    except OSError as e:
        print(e)
        pass

    # else, try to copy it to current folder and run it directly
    shutil.copy2(blast_path+'rpsblast','rpsblast')
    rpsblast_cline = NcbirpsblastCommandline('rpsblast', query=blast_file, db=rpsblast_db, evalue=0.00000001,outfmt=5, out="rpsblastOutput.xml")
    return rpsblast_cline() 

def getSubjectInfo(blast_records, prots_of_interest, evalue):
    '''Run through an iterator of blast records and save the results to a dictionary
    
       Args:
           blast_records (iterable):	An iterable of blast records
           prots_of_interest (dict):	A dictionary with proteins of interest. 
                                                                            This contains as key protein name (has to be same as protein names found in blast_records)
                                                                            and as value inform
    
       Returns:
           Dictionary with as key name of the subject (found) protein of a BLAST records and as value a dict with as key query (searched) protein from BLAST records
    '''
    subject_info = {}
    
    #see if target is in TF list
    for blast_record in blast_records:
        query = blast_record.query
        for alignment in blast_record.alignments:
            subj_title = alignment.title
            # for some reason makeblastdb can put BL_ORD_id in front of the actual description
            if re.match('\Agnl\|BL_ORD_ID\|\d+\s',subj_title) > 0:
                subj_title = re.split('\Agnl\|BL_ORD_ID\|\d+\s',subj_title)[1]
            isInterest = False
            for hsp in alignment.hsps:
                # check if the protein is not a protein of interest (the ones used for query) and evalue of blast lower than given evalue
#                if 'zpr' in subj_title.lower():
#                    print subj_title, hsp.expect,evalue, float(hsp.expect) <= float(evalue)
                if float(hsp.expect) <= float(evalue):
#                    if 'zpr' in subj_title.lower():
#                        print 'got added to subject_info'
                    #saving all unique genes and their info from tair db
                    if subject_info.has_key(subj_title):
                        subject_info[subj_title]['hsp'].append(hsp)
                        subject_info[subj_title]['query_title'].append(query)
                    else:
                        subject_info[subj_title] = {'hsp':[hsp], 'query_title':[query], 'domains':[]}
    return subject_info

def subjectHitsInterest(blast_records, proteins_of_interest, evalue):
    '''Find the subject (used to search) proteins in blast_records that are also in proteins_of_interest
    
       Args:
               subject_info (dict):			   Key is name of subject protein from previous BLAST, value is info query protein from previous BLAST
               proteins_of_interest (dict):	   Key is name of proteins of interest (e.g. transcription factors)
               blast_records (dict):		   Iterator of blast records
               evalue (float):				   Minimum e-value to be included
       Returns:
               A dict like subject_info, but filtered for proteins that found a protein of interest in blast_records
    '''
    # found_interest will contain all the proteins from subject_info that had a hit with a protein of interest
    subject_info = {}
    for blast_record in blast_records:
#        to_print = False
#        if 'zpr' in blast_record.query.lower():
#            to_print = True
#            print blast_record.query
        for alignment in blast_record.alignments:
            # BLAST adds an gnl|BL_ORD_ID|<number> to the titles. needs to be removed
            if re.match('\Agnl\|BL_ORD_ID\|\d+\s',alignment.title) > 0:
                alignment.title = re.split('\Agnl\|BL_ORD_ID\|\d+\s',alignment.title)[1]
        
            for hsp in alignment.hsps:
#                if to_print:
#                    print alignment.title
#                    print hsp.expect, evalue, proteins_of_interest.has_key(str(alignment.title)), float(hsp.expect) <= float(evalue)
                # aligment title is the title of the hits from the subject_info (subjects from subject_info are queries in this BLAST search)
                if proteins_of_interest.has_key(str(alignment.title)) and float(hsp.expect) <= float(evalue):
                    if subject_info.has_key(blast_record.query):
                        subject_info[blast_record.query]['hsp'].append(hsp)
                        subject_info[blast_record.query]['query_title'].append(alignment.title)
                    else:
                        subject_info[blast_record.query] = {'hsp':[hsp], 'query_title':[alignment.title], 'domains':[]}
                            
    return subject_info

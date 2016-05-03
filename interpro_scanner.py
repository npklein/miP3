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

import iprscan_soappy
import re
import pickle
import os.path
import traceback
import SOAPpy
import time

def interpro_result(interpro_submit_sequences, email, developing):
    protein_ipr_db_domain = {}
    # this is done per 25
    for protein_name, interpro_result in iprscan_soappy.runInterpro(interpro_submit_sequences, email):					# get dict with as value protein name and as value various stuff
        ipr_domain_names = []
        protein_ipr_db_domain[protein_name] = {}

        for ipr_code in interpro_result:
            # list of ipr domain names for this protein
            if 'ipr_names' in interpro_result[ipr_code]:
                ipr_domain_names += interpro_result[ipr_code]['ipr_names']
            for database in interpro_result[ipr_code]:
                protein_ipr_db_domain[protein_name][ipr_code] = {database:interpro_result[ipr_code][database]}	 # update it with database and database specific name
                # make a separate list for PFAM domains, because these are used later
        protein_ipr_db_domain[protein_name]['ipr_domain_names'] = ipr_domain_names
        if developing:
            try:
                interpro_data = os.path.dirname(os.path.abspath(__file__))+os.sep+'interpro_results'+os.sep+protein_name.split('|')[0].strip()+'_interpro.p'
                f = open( interpro_data, 'wb' )
                pickle.dump( protein_ipr_db_domain[protein_name], f )
                print 'wrote interpro data to '+interpro_data
            except:
                print traceback.print_exc()
                pass
    return protein_ipr_db_domain

def interproScan(subject_info, all_proteins,pfam_domains_file, email, developing):
    '''Scans all the proteins in protein_list with InterproScan. Remove subjects from subject_info that have a domain that is now allowed (as given by user, e.g. transcription factor domains)
        or that don't have any IPR domains in comon with any of its query proteins.
        
        Args:
        subject_info (dict):	   Dict with as key name of subject (found) protein and as one of the values the name of the query (searched) protein(s)
        all_proteins (dict):	   Dict that contains as key protein name and as value the sequence of that protein
        email (string):            email adres that will be send to interproscan
        
        Returns:
        Dict that is subject_info but filtered on subjects from subject_info that have a domain that is now allowed (as given by user, e.g. transcription factor domains)
        or that don't have any IPR domains in comon with any of its query proteins.
        '''
    sequence_count = 0														  # sequence counter for status bar update
    fasta_sequences = set()													 # use set so no double sequence are scanned
    # get all the query proteins
    for sequence in subject_info:
        fasta_sequences.add('>'+sequence+'\n'+str(all_proteins[sequence]))			 # add fasta format string of all subject (found) proteins
        for query in subject_info[sequence]['query_title']:
            fasta_sequences.add('>'+query+'\n'+str(all_proteins[query]))			 # add fasta format string of all query (searched) proteins

    interpro_submit_sequences = []	# submit 15 asynchronyous jobs from this list
    length_of_sequence_list = len(fasta_sequences)								 # needed to update the statusbar a little bit per sequence
    if length_of_sequence_list == 0:	  # the length added per sequence. 
        print('Nothing found','No putative mip found, exiting program')								 # tell about the error
    x=1
    protein_ipr_db_domain = {}
    for fasta_sequence in fasta_sequences:											 # loop over all the sequences
        print(str(x)+'/'+str(len(fasta_sequences)))
        x+=1
        sequence_count += 1
        protein_name = (fasta_sequence.split('\n')[0].lstrip('>'))					 # get the protein name
        if not os.path.exists('interpro_results'):
            os.makedirs('interpro_results')
        interpro_file = os.path.dirname(os.path.abspath(__file__))+os.sep+'interpro_results'+os.sep+protein_name.split('|')[0].strip()+'_interpro.p'
        if developing:
            if os.path.isfile(interpro_file):
                try:
                    interpro_data = pickle.load( open(interpro_file, 'rb' ) )
                    protein_ipr_db_domain[protein_name] = interpro_data
                    print 'loaded interpro data from '+str(interpro_file)
                except EOFError as e:
                    print(e)
                    print('removing pickle file, blasting again')
                    # if file is not read correctly, remove it so that it is remade
                    os.remove(interpro_file)
        if not os.path.isfile(interpro_file):
            # simple counter
            protein_ipr_db_domain[protein_name] = {}									# keep the ipr codes and database name and domain name in this dict
            interpro_submit_sequences.append(fasta_sequence)
            if sequence_count % 25 == 0 or len(fasta_sequences) <= sequence_count-1:	  # divivide it up into chunks of 25
                y = 0
                while True:
                    try:
                        interpro_data = interpro_result(interpro_submit_sequences,email,developing)
                        break
                    except SOAPpy.Errors.HTTPError:
                        y+=1
                        if y == 100:
                            print('failed 100 times, something wrong with interpro server')
                            raise
                        print 'HTTPError, sleep 1 minute to give running jobs time to finish, then submitting same jobs again'
                        print 'try '+str(y)+'/100'
                        time.sleep(60)
                        
                protein_ipr_db_domain.update(interpro_data)
                interpro_submit_sequences = []

    # if for some reason not all of them were done
    if len(interpro_submit_sequences) > 0:
        interpro_data = interpro_result(interpro_submit_sequences, email, developing)
        protein_ipr_db_domain.update(interpro_data)

#### remove subjects with a pfam domain or IPR code that's not allowed and subjects that don't have any domains in common with their query protein(s)
    subject_info_filtered_on_domains = {}											# First, remove the proteins that have a Pfam domain that's not allowed and sve in this dict
    pfam_domains_from_user_file = open(pfam_domains_file).read()  # Read the pfam domains file as given by user
    unwanted_domain_set = set(re.split('\n',pfam_domains_from_user_file)) # create a set by splitting (splits on enters;
    unwanted_domain_set = [x.lower().strip() for x in unwanted_domain_set]			# lower the domains and strip all whitespaces, so they can be matched with other domains
#    print 'the following domains get filtered out:'
#    print unwanted_domain_set
    for subject in subject_info:
        # filter on unwanted pfam domains again, even if it was done before, because server is most up to date. Previous pfam filter was to
        # pre-eliminate proteins so less have to go through interproscan.
        if subject not in protein_ipr_db_domain:
#            try:
            print('missed getting interpro file for',subject, ', searching interpro now')
            x = 0
            while True:
                try:
                    result = interpro_result(['>'+subject+'\n'+str(all_proteins[subject])])
                except SOAPpy.Errors.HTTPError:
                    x+=1
                    if x == 100:
                        print('failed 100 times, something wrong with interpro server')
                        raise
                    print 'HTTPError, sleep 1 minute to give running jobs time to finish, then submitting same jobs again'
                    time.sleep(60)
            protein_ipr_db_domain.update(result)
#            except:
#                print subject, 'could not be processed'
#                continue
        ipr_names = [i for i in protein_ipr_db_domain[subject].keys() if i != 'ipr_domain_names']
        domains = protein_ipr_db_domain[subject]['ipr_domain_names'] + ipr_names
#        if 'AT1G72416.3' in subject:
#            print 'SUBJECT:',subject
#            print 'protein_ipr_db_domain',set(protein_ipr_db_domain[subject]['ipr_domain_names'])
#            print domains
        if set(domains).isdisjoint(unwanted_domain_set):# and len(protein_domains[subject]) > 0:				   # if domains from subject are not in domain_set: check if they have at leat 1 domain the query also has																		 # if none of the domains are in the unwanted pfam file
            # remove the proteins that don't have IPR ids in common with their 'query' proteins, or if query protein doesn't have any domains
            #            print set(protein_domains[subject])
            
            for query in subject_info[subject]['query_title']:
#                if 'AT1G72416.3' in subject:
#                    print 'SUBJECT:',subject
#                    print 'QUERY:',query
#                    print 'protein_ipr_db_domain[subject][\'ipr_domain_names\']:', protein_ipr_db_domain[subject]['ipr_domain_names']
#                    print 'protein_ipr_db_domain[query][\'ipr_domain_names\']:',protein_ipr_db_domain[query]['ipr_domain_names']
#                    print 'not set(protein_ipr_db_domain[subject][\'ipr_domain_names\']).isdisjoint(protein_ipr_db_domain[query][\'ipr_domain_names\']):',not set(protein_ipr_db_domain[subject]['ipr_domain_names']).isdisjoint(protein_ipr_db_domain[query]['ipr_domain_names'])
                if not set(protein_ipr_db_domain[subject]['ipr_domain_names']).isdisjoint(protein_ipr_db_domain[query]['ipr_domain_names']):	  # if there is overlap between the subject IPR IDs and the query IPR domains
                    subject_info_filtered_on_domains[subject] = subject_info[subject]	 #	  add to the dict
                    subject_info_filtered_on_domains[subject]['domains'] = protein_ipr_db_domain[subject]['ipr_domain_names']
                elif len(protein_ipr_db_domain[subject]['ipr_domain_names']) == 0: # add proteins without domains
                    subject_info_filtered_on_domains[subject] = subject_info[subject]	 #	  add to the dict
                    subject_info_filtered_on_domains[subject]['domains'] = ''
            if subject in subject_info_filtered_on_domains:
                subject_info_filtered_on_domains[subject]['ipr_domains'] = ipr_names
    del(subject_info)															# delete for memory
    print('len subject_info_total after removing proteins with unwanted pfam domains \nand which have none in common with their query protein(s): '+str(len(subject_info_filtered_on_domains))+'\n')
    return subject_info_filtered_on_domains


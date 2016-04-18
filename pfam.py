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

def pfamSearch(subject_info, all_proteins,pfam_domain_file):
    '''Searches against Pfam database. miP.py uses InterPro scan instead as main tool, but this is used for speed reasons. 
    '''
    pfam_search_str = ''
    protein_set = set()
    for subject in subject_info:
        if subject not in protein_set:
            pfam_search_str += '>'+subject+'\n'+str(all_proteins[subject])+'\n'
            protein_set.add(subject)
    pfam_search_file = open('prots_to_get_domains.fasta','w')
    pfam_search_file.write(pfam_search_str)
    pfam_search_file.close()
    blast_records = blast.rpsblast(self.path+'prots_to_get_domains.fasta', self.path+'Pfam', self.path, evalue=1)
    # filter out subjects that have the domain to be filtered out
    protDomains = domains.getDomains(blast_records)
    pfam_domains_from_user_file = open(self.pfam_domains_file).read()  # Read the pfam domains file as given by user
    unwanted_domain_set = set(re.split(';|,|\t|\n|\s',pfam_domains_from_user_file)) # create a set by splitting (splits on comma, whitespace, enters, tabs and ;
    unwanted_domain_set = [x.lower().strip() for x in unwanted_domain_set] 
    subject_info_filter_on_unwanted_domains = {}
    for subject in protDomains:
        if protDomains[subject].isdisjoint(unwanted_domain_set): 
            subject_info_filter_on_unwanted_domains[subject] = subject_info[subject]
    # add proteins that did not have any domains
    for subject in subject_info:
        if not protDomains.has_key(subject):
            subject_info_filter_on_unwanted_domains[subject] = subject_info[subject]
    sys.stderr.write('len subject_info after filtering on unwanted PFAM domains: '+str(len(subject_info_filter_on_unwanted_domains))+'\n')
    return subject_info_filter_on_unwanted_domains

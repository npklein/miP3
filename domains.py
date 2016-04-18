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

import blast
import os


def getDomains(blast_records):   
    '''Get the domains out of Pfam BLAST results
    
       Args:
           blast_records (an iterator of a Blast record):    BLAST records of a PFAM BLAST search
        
       Returns:
           Dict[query protein] = [list, of, domains].
    '''
    protDomains = {}
    for blast_record in blast_records:
    #We want to ignore any queries with no search results:
        if blast_record.alignments:
            for alignment in blast_record.alignments :
                subjectDomain = alignment.title.split(".")[0].split(",")[1].replace(" ","").lower()
                if protDomains.has_key(blast_record.query):
                    protDomains[blast_record.query].add(subjectDomain)
                else:
                    protDomains[blast_record.query] = set([subjectDomain])
    return protDomains
    

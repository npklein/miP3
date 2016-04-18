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

def write_fasta(protein_ids, all_proteins,  outfile_location):
    subject_str = ''
    for protein in protein_ids:
        try:
            subject_str += '>'+protein+'\n'+str(all_proteins[protein])+'\n'
        except KeyError:
            print 'Error:',protein,'not in the all protein file.'
            exit(-1)
    out_file = open(outfile_location,'w')
    out_file.write(subject_str)
    out_file.close()
def write_result_out(subject_info, outfile_location, all_proteins):
    out_str = 'miP\tTFs\tmiP domains\tmiP length\n'
    written = []
    for mip in subject_info:
        if not mip in written:
            out_str += mip+'\t'+', '.join(subject_info[mip]['query_title'])+'\t'+', '.join(subject_info[mip]['ipr_domains'])+'\t'+str(len(all_proteins[mip]))+'\n'
        written.append(mip)
    out_file = open(outfile_location,'w')
    out_file.write(out_str)
    out_file.close()
    print('result written to: '+str(outfile_location))

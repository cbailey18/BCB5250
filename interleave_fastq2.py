# -*- coding: utf-8 -*-
#!/usr/bin/env python3
import argparse
#from Bio import SeqIO
from string import Template
import os

# Definte the function interleave_fastq, which will take in 2 fastq files and output an interleaved fastq file

def interleave_fastq(forward_file, reverse_file, output_file):

# Open the forward and reverse files
    fwd_file = open(forward_file)
    rev_file = open(reverse_file)

# Read the files and store each line as an element in a list
    fwd_read = fwd_file.readlines()
    rev_read = rev_file.readlines()

# Strip the /n characters added to each line during the readline process
    for i in range(0, len(fwd_read)):
        fwd_read[i] = fwd_read[i].strip()

    for i in range(0, len(rev_read)):
        rev_read[i] = rev_read[i].strip()


  # create a list of lists where each list contains lines x to x+4 from the input files
    fwd_byread = [fwd_read[x:x+4] for x in range(0, len(fwd_read)) if len(fwd_read[x:x+4])]
    rev_byread = [rev_read[x:x+4] for x in range(0, len(rev_read)) if len(rev_read[x:x+4])]

  # filter out any sublist that does not start with '@' to narrow files to only those starting with sequence id
    fwd_output = [sublist for sublist in fwd_byread if sublist[0][0] =='@']
    rev_output = [sublist for sublist in rev_byread if sublist[0][0] =='@']

# create an empty list and append fwd and rev reads at the same index in alternating order
 # interleave_list = []
  #for i in range(0, len(fwd_output)):
   # interleave_list.append(fwd_output[i])
    #interleave_list.append(rev_output[i])

  # create a template for the fastq format
    fastq_template = Template(f'$identifier\n$sequence\n$plus_sign\n$quality_score\n')

    # Ensure directory for output exists before writing to it
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

  # open the output file so that it can be written to. Iterate through the previously created interleaved list. Take each element of each sublist and add it to its own line in the output file
    with open(file = output_file, mode='w') as file:
        for i in range(0, len(fwd_output)):
            identifier_fwd = fwd_output[i][0]
            sequence_fwd = fwd_output[i][1]
            sign_fwd = fwd_output[i][2]
            quality_score_fwd = fwd_output[i][3]
            file.writelines(fastq_template.substitute(identifier=identifier_fwd,
                                         sequence=sequence_fwd,
                                         plus_sign=sign_fwd,
                                         quality_score=quality_score_fwd))
            identifier_rev = rev_output[i][0]
            sequence_rev = rev_output[i][1]
            sign_rev = rev_output[i][2]
            quality_score_rev = rev_output[i][3]
            file.writelines(fastq_template.substitute(identifier=identifier_rev,
                                         sequence=sequence_rev,
                                         plus_sign=sign_rev,
                                         quality_score=quality_score_rev))

  #close the input files
    fwd_file.close()
    rev_file.close()

    if os.path.exists(output_file):
        print("yes")
    else:
        print("no")
  


# This code included in provided template but not used
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Interleave Paired-End FASTQ Files')
    parser.add_argument('forward_file', help = 'Forward reads FASTQ file')
    parser.add_argument('reverse_file', help = 'Reverse reads FASTQ file')
    parser.add_argument('out_file', help = 'Output interleaved FASTQ file')
    args = parser.parse_args()
    interleave_fastq(args.forward_file, args.reverse_file, args.out_file)



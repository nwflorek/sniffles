#!/usr/bin/env python3
#Author Kelsey Florek

import re
import os,sys

path = sys.argv[1]

fasta_files = []
vcf_files = []
for root,dirs,files in os.walk(path):
    for file in files:
        if '.fasta' in file:
            fasta_files.append(file)
        if '.vcf' in file:
            vcf_files.append(file)

for fasta_file in fasta_files:
    fasta_seq = SeqIO.parse(open(os.path.join(path,fasta_file)),'fasta')
    for fasta in fasta_seq:
        outname = fasta_file.split('.')[0]
        fasta.id = outname
        with open(os.path.join(outdir,outname+'.fasta'),'w') as output_fasta:
            SeqIO.write(fasta,output_fasta,'fasta')

for vcf_file in vcf_files:
    lines = []
    id = vcf_file.split('_')[0]
    with open(os.path.join(path,vcf_file),'r') as infile:
        lines = infile.readlines()

    current_name = ''
    for line in lines:
        if line[0] != '#':
            current_name = line.split('\t')[0]
            break

    with open(os.path.join(path,vcf_file),'w') as outfile:
        for line in lines:
            new_line = re.sub(re.escape(current_name),re.escape(id), line)
            outfile.write(new_line)

#!/usr/bin/env python3
#Author Kelsey Florek
#Description: This script takes the freqs generated from freebayes
#(freebayes -f ref.fa --haplotype-length 0 --min-alternate-count --min-alternate-fraction 0 --pooled-continuous --report-monomorphic >var.vcf)
#and generates a allele matrix that can be used in downstream analyses.

import re
import os,sys
from cyvcf2 import VCF
import csv

path = sys.argv[1]
vcf = VCF(path)
samples = vcf.samples

class _var():
    def __init__(self,sample,chrom,pos,ref,depth):
        self.sample = sample
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.depth = depth
        self.freqs = {'A':0,'T':0,'G':0,'C':0}
    def set_freq(self,base,freq):
        if self.freqs[base] != 0.0:
            print(self.sample + " Already has a a freq of " + str(self.freqs[base]) + " at " + self.chrom + ' ' + str(self.pos))
        else:
            self.freqs[base] = freq
    def print_var(self):
        print(self.sample,self.chrom,self.pos,self.ref,self.depth,self.freqs)

variants = {}
for sample in samples:
    vcf = VCF(path)
    vcf.set_samples(sample)
    for var in vcf:
        dp = var.format("DP")[0]
        ad = var.format("AD")[0]
        alts = var.ALT
        bases = [var.REF]
        for alt in alts:
            bases.append(alt)
        name = sample +str(var.CHROM)+str(var.POS)
        if dp[0] != -2147483648:
            i = 0
            while i < len(bases):
                if name not in variants:
                    variants[name] = _var(sample,var.CHROM,var.POS,var.REF,dp[0])
                    variants[name].set_freq(bases[i],ad[i])
                else:
                    variants[name].set_freq(bases[i],ad[i])
                i += 1

for key in variants:
    variants[key].print_var()

with open('variant_list.csv','w',newline='') as csvfile:
    csvwriter = csv.writer(csvfile,delimiter=',')
    for key in variants:
        vars = [variants[key].freqs['A'],variants[key].freqs['T'],variants[key].freqs['G'],variants[key].freqs['C']]
        csvwriter.writerow([variants[key].sample,variants[key].chrom+"_"+str(variants[key].pos),vars])

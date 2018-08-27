#!/usr/bin/env python3
# Sniffles
#Author: Nick Florek
#Performs SNP analysis of influenza genomes from raw reads.
import yaml
import argparse
import os
from shutil import copyfile
import multiprocessing as mp
from trim import trimmomatic
from mapping import mapping
from consensus import consensus
import readcleaning as rc
from snpcaller import snpcaller
import filehandler as fh

#determine command line arguments and get path
parser = argparse.ArgumentParser(description='Pipeline to examine SNPs from raw illumina reads')
parser.add_argument('-c',metavar='config', type=str,help="config file")
parser.add_argument('-i',metavar='input', type=str,help="raw reads directory - defaults to working directory")
parser.add_argument('-o',metavar='output', type=str,help="output directory - defaults to working directory")
parser.add_argument('-t',metavar='threads', type=int,help="number of cpus to use for pipeline",default=4)
args = parser.parse_args()
numThreads = args.t
configFile = args.c
libPath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'lib')

#get input and output paths
try:
    outDir = os.path.abspath(args.o)
except (AttributeError, TypeError) as err:
    outDir = os.getcwd()
try:
    inDir = os.path.abspath(args.i)
except (AttributeError, TypeError) as err:
    inDir = os.getcwd()

#open config file and store configuation
with open(configFile,'r') as ymlFile:
    cfg = yaml.load(ymlFile)

#create outdir
try:
    os.mkdir(os.path.join(outDir,cfg['exec']['outdir']))
except FileExistsError:
    print('Please remove the files generated from the previous run or rename the output file.')
    exit()

#copy reference to outdir
copyfile(cfg['exec']['referenceSequence'],os.path.join(outDir,cfg['exec']['outdir'])+'/'+cfg['exec']['referenceSequence'])

#parse and store read information from input directory
r = fh.reads(inDir)

#begin trimming
pool = mp.Pool(processes=numThreads)
pool.starmap(trimmomatic,[[libPath,cfg,i[0],i[1],i[2]] for i in r.retList()])

#begin mapping
mapping(libPath,cfg,numThreads,r.idList)

#datacleaning
if cfg['exec']['removeDupReads']:
    rc.removeDuplicates(libPath,cfg,numThreads,r.idList)

#normalize coverage
if cfg['exec']['normalizeCoverage']:
    rc.normCoverage(libPath,cfg,numThreads,r.idList)

#generate consensus
if cfg['exec']['generateConsensus']:
    consensus(libPath,cfg,numThreads,r.idList)

#call snps
if cfg['exec']['callSNPs']:
    snpcaller(libPath,cfg,numThreads,r.idList)

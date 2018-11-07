#!/usr/bin/env python3
# Sniffles
#Author: Nick Florek
#Performs SNP analysis of influenza genomes from raw reads.
import yaml
import argparse
import os
from shutil import copyfile
from trim import trimmomatic
from mapping import mapping
from consensus import consensus
import readcleaning as rc
from snpcaller import snpcaller
import filehandler as fh
import time
import datetime
import sc

#determine command line arguments and get path
parser = argparse.ArgumentParser(description='Pipeline to examine SNPs from raw illumina reads')
parser.add_argument('-c',metavar='config', type=str,help="config file")
parser.add_argument('-i',metavar='input', type=str,help="raw reads directory - defaults to working directory")
parser.add_argument('-o',metavar='output', type=str,help="output directory - defaults to working directory")
parser.add_argument('-t',metavar='threads', type=int,help="number of cpus to use for pipeline",default=4)
args = parser.parse_args()
numThreads = args.t
configFile = args.c

#print main display title
sc.mainTitle()

#get start time
start = time.time()

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

#add lib path to cfg
cfg['libPath'] = os.path.join(os.path.dirname(os.path.abspath(__file__)),'lib')

#create outdir
cfg['exec']['outdir'] = os.path.join(outDir,cfg['exec']['outdir'])
try:
    os.mkdir(cfg['exec']['outdir'])
except FileExistsError:
    print('Please remove the files generated from the previous run or rename the output file.')
    exit()

#copy reference to outdir
copyfile(cfg['exec']['referenceSequence'],os.path.join(outDir,cfg['exec']['outdir'])+'/'+cfg['exec']['referenceSequence'])

#parse and store read information from input directory
readData = fh.reads(inDir)

trimmomatic(readData,cfg,numThreads)

#begin mapping
mapping(readData,cfg,numThreads,jobtype='map-trimmed')
#report maping average depth
with open(os.path.join(cfg['exec']['outdir'],'average_depth.log'),'w') as outdepth:
    for key in readData.data['mapData']['avgDepth']:
        avgdepth = readData.data['mapData']['avgDepth'][key]
        outdepth.write(f'{key},{avgdepth}\n')

#datacleaning
#normalize coverage
if cfg['exec']['normalizeCoverage']:
    rc.normCoverage(readData,cfg,numThreads)

#remove duplicate reads
if cfg['exec']['removeDupReads']:
    rc.removeDuplicates(readData,cfg,numThreads)

#generate consensus
if cfg['exec']['generateConsensus']:
    consensus(readData,cfg,numThreads)

#call snps
if cfg['exec']['callSNPs']:
    snpcaller(readData,cfg,numThreads)

sc.procTitle('Finished Sniffles')
end = time.time()
runtime = round(end - start,2)
runtime = str(datetime.timedelta(seconds=runtime))
print(f'Sniffles finished with a total runtime of {runtime}.')

#!/usr/bin/env python3
# Sniffles
#Author: Kelsey Florek
#Performs SNP analysis of influenza genomes from raw reads.
import yaml
import argparse
import os, sys
from shutil import copyfile
from trim import trim
from mapping import mapping,indexing,average_depth
from consensus import consensus
import readcleaning as rc
from snpcaller import snpcaller
import fileparser as fp
import time
import datetime
import sc

#print main display title
sc.mainTitle()

#determine command line arguments and get path
parser = argparse.ArgumentParser(description='Pipeline to examine SNPs from raw illumina reads')
parser.add_argument('-c',metavar='config', type=str,help="config file",default='config.yml')
parser.add_argument('-i',metavar='input', type=str,help="raw reads directory - defaults to working directory")
parser.add_argument('-o',metavar='output', type=str,help="output directory - defaults to working directory")
parser.add_argument('-t',metavar='threads', type=int,help="number of cpus to use for pipeline",default=4)
if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
args = parser.parse_args()
numThreads = args.t
configFile = args.c

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

#create outdir
cfg['exec']['outdir'] = os.path.join(outDir,cfg['exec']['outdir'])
try:
    os.mkdir(cfg['exec']['outdir'])
except FileExistsError:
    cfg['exec']['outdir'] = cfg['exec']['outdir']+'_'+str(int(time.time()))
    os.mkdir(cfg['exec']['outdir'])


#copy reference to outdir
copyfile(cfg['exec']['referenceSequence'],os.path.join(outDir,cfg['exec']['outdir'])+'/'+cfg['exec']['referenceSequence'])

#parse and store read information from input directory
readData = fp.RunFiles(inDir)

#trim the reads
trim(readData,cfg,numThreads)

#setup inital mapping jobs
mapping_list = []
for id in readData.runtime['trimmed']:
    mapping_list.append((id,readData.runtime['trimmed'][id][0],readData.runtime['trimmed'][id][1],os.path.abspath(cfg['exec']['referenceSequence'])))
sc.checkexists(os.path.join(cfg['exec']['outdir']+'/inital_mapping'))


#index reference reference sequence
indexing(cfg,os.path.abspath(cfg['exec']['referenceSequence']))

#run inital mapping jobs
bam_list = mapping(cfg,mapping_list,cfg['exec']['outdir']+'/inital_mapping',numThreads)

#determine average depth
print('\nSniffles: Calculating average read depth in the initial mapping')
with open(os.path.join(cfg['exec']['outdir'],'average_depth.log'),'w') as outdepth:
    filtered_bam_list = []
    for path in bam_list:
        filename = os.path.basename(path)
        id = filename.split('.')[0]
        avgdepth = average_depth(path)

        #only add isolates that pass average depth
        if float(avgdepth) >= cfg['exec']['minimumAverageDepth']:
            outdepth.write(f'{id},{avgdepth} : Pass\n')
            filtered_bam_list.append(path)
        else:
            outdepth.write(f'{id},{avgdepth} : Fail\n')
    bam_list = filtered_bam_list
print('\nSniffles: Finished determining average read depth')

#datacleaning

#normalize coverage
if cfg['exec']['normalizeCoverage']:
    fastq_list = rc.normCoverage(cfg,bam_list,numThreads)
    mapping_list = []
    for fastq in fastq_list:
        id = os.path.basename(fastq).split('.')[0]
        mapping_list.append((id,fastq,'',os.path.abspath(cfg['exec']['referenceSequence'])))
    sc.checkexists(os.path.join(cfg['exec']['outdir']+'/norm_mapping'))
    bam_list = mapping(cfg,mapping_list,cfg['exec']['outdir']+'/norm_mapping',numThreads)

#remove duplicate reads
if cfg['exec']['removeDupReads']:
    bam_list = rc.removeDuplicates(cfg,bam_list,numThreads)

#generate consensus
if cfg['exec']['generateConsensus']:
    fasta_list = consensus(cfg,bam_list,numThreads)

    #map reads to consensus
    if cfg['exec']['mapToConsensus']:
        mapping_list = []
        indexing(cfg,*fasta_list)
        for id in readData.runtime['trimmed']:
            for fasta in fasta_list:
                fasta_id = os.path.basename(fasta).split('.')[0]
                if fasta_id == id:
                    mapping_list.append((id,readData.runtime['trimmed'][id][0],readData.runtime['trimmed'][id][1],os.path.abspath(fasta)))
        mapping(cfg,mapping_list,cfg['exec']['outdir']+'/map_to_consensus',numThreads)


#call snps
if cfg['exec']['callSNPs']:
    snpcaller(cfg,bam_list,numThreads)

sc.procTitle('Finished Sniffles')
end = time.time()
runtime = round(end - start,2)
runtime = str(datetime.timedelta(seconds=runtime))
print(f'Sniffles: Finished with a total runtime of {runtime}.')

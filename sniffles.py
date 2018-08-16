#!/usr/bin/env python3
# Sniffles
#Author: Nick Florek
#Performs SNP analysis of influenza genomes from raw reads.
import yaml
import argparse
import os
import multiprocessing as mp
from trim import trimmomatic
from mapping import mapping

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

#class for containing read information
class reads:
    #list contianing all isolate ids
    idList = []
    #data dictonary containing all of the run information
    #formatted id: ["read1 path","read2 path"]
    data = {}
    def __init__(self,path):
        #list of reads
        readList = []
        for root,dirs,files in os.walk(path):
            #scan path and look for fastq files and record ids
            for file in files:
                if '.fastq' in file:
                    if '_R1' in file or '_1' in file:
                        id = file.split('_')[0]
                        if id not in self.idList:
                            self.idList.append(id)
                    if '_R2' in file or '_2' in file:
                        id = file.split('_')[0]
                        if id not in self.idList:
                            self.idList.append(id)
                    readList.append(root+'/'+file)
        readList.sort()
        for id in self.idList:
            for read in readList:
                if id in read and '_R1' in read or '_1' in read:
                    self.data[id] = [read]
            for read in readList:
                if id in read and '_R2' in read or '_2' in read:
                    self.data[id].append(read)

    #return a list of id with each paired path as a 3 item sublist
    def retList(self,):
        l = []
        for id in self.idList:
            l.append([id,self.data[id][0],self.data[id][1]])
        return l

r = reads(inDir)

#begin trimming
pool = mp.Pool(processes=numThreads)
pool.starmap(trimmomatic,[[libPath,cfg,i[0],i[1],i[2]] for i in r.retList()])

#begin mapping
mapping(libPath,cfg,numThreads,r.idList)

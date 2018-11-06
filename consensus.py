import os
import shlex
import subprocess as sub
import time
from sniffProc import proc,init
from sc import procTitle,checkexists
from mapping import mapping
import multiprocessing as mp

def consensus(readData,runCFG,threads='1',ids=''):
    #easy const access
    libPath = runCFG['libPath']
    outDir =runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])

    #if no id supplied get list of all ids
    if not ids:
        ids = readData.idList

    #if no consensus dir exists create one
    try:
        os.mkdir(os.path.join(outDir,'consensus'))
    except:
        pass

    #notify starting mapping
    procTitle('Generate Consensus')
    print('\nSniffles: Started generating consensus sequence')
    #get start time
    overall_start = time.time()
    start = time.time()
    #set reference sequence and begin
    reference = runCFG['exec']['referenceSequence']

    #setup multiprocessing
    lock = mp.Lock()
    pool = mp.Pool(processes=threads,initializer=init,initargs=(lock,))

    #command list for generating mpileups
    cmds = []
    sort_cmds = []
    for id in ids:
        #determine samfile that will be used
        #no read cleaning
        if runCFG['exec']['removeDupReads']:
            samfile = f'{outDir}/nodups/{id}.sam'
        elif runCFG['exec']['normalizeCoverage']:
            samfile = f'{outDir}/normalized_mapping/{id}.sam'
        else:
            samfile = f'{outDir}/inital_mapping/{id}.sam'

        #convert sam to sorted bam file
        sort_cmd = f'{libPath}/bin/samtools view -b {samfile} | samtools sort > {outDir}/consensus/{id}.bam'
        sort_cmds.append(sort_cmd)

        #make multiway pileup using samtools
        cmd = f'{libPath}/bin/samtools mpileup -d 1000000 {outDir}/consensus/{id}.bam -f {reference} -o {outDir}/consensus/{id}.pileup'
        cmds.append(cmd)

    #start multiprocessing
    pool.starmap(proc, [[runCFG,i,'','',True] for i in sort_cmds])
    pool.starmap(proc, [[runCFG,i] for i in cmds])

    #remove temp bam files
    for id in ids:
        os.remove(f'{outDir}/consensus/{id}.bam')

    #remove pileup file if it is empty
    pile_ids = []
    for id in ids:
        try:
            file = os.path.join(outDir,'consensus',f'{id}.pileup')
            if not os.path.getsize(file)>0:
                os.remove(file)
            else:
                pile_ids.append(id)
        except:
            pass
    ids = pile_ids

    end = time.time()
    runtime = round(end - start,2)
    print(f'\nSniffles finished pileup in {runtime} seconds')
    start = time.time()

    #command list for generating consensus vcf
    cmds = []
    outFiles = []
    for id in ids:
        #run varscan mpileup2cns to generate vcf with consensus information
        minCov = runCFG['snpcalling']['minCoverage']
        quality = runCFG['snpcalling']['snpQualityThreshold']
        freq = runCFG['snpcalling']['consensusFrequency']
        cmd = f'java -jar {libPath}/varscan/VarScan.v2.3.9.jar mpileup2cns {outDir}/consensus/{id}.pileup --min-coverage {minCov} --min-avg-qual {quality} --min-var-freq {freq} --strand-filter 1 --output-vcf 1'
        cmds.append(cmd)
        o = f'{id}.vcf'
        outFiles.append(o)

    #start multiprocessing
    pool.starmap(proc, [[runCFG,cmds[i],'consensus',outFiles[i]] for i in range(len(cmds))])

    #check if vcf file is empty, if it is skip id and remove vcf file
    vcf_ids = []
    for id in ids:
        try:
            vcf_file = os.path.join(outDir,'consensus',f'{id}.vcf')
            if os.path.getsize(vcf_file)>0:
                empty = True
                with open(vcf_file,'r') as vcffile:
                    for line in vcffile:
                        line = line.strip()
                        if line and line[0] != '#':
                            empty = False
                if empty:
                    os.remove(vcf_file)
                else:
                    vcf_ids.append(id)
            else:
                os.remove(vcf_file)
        except:
            pass
    ids = vcf_ids

    end = time.time()
    runtime = round(end - start,2)
    print(f'\nSniffles finished generating the VCF in {runtime} seconds')
    start = time.time()

    #command list for compressing files
    zip_cmds = []
    idx_cmds = []
    for id in ids:
        #compress vcf file with bgzip
        cmd = f'{libPath}/bin/bgzip {id}.vcf'
        zip_cmds.append(cmd)
        #index compressed vcf with tabix
        cmd = f'{libPath}/bin/tabix {id}.vcf.gz'
        idx_cmds.append(cmd)
    #start multiprocessing
    pool.starmap(proc, [[runCFG,i,'consensus'] for i in zip_cmds])
    pool.starmap(proc, [[runCFG,i,'consensus'] for i in idx_cmds])

    end = time.time()
    runtime = round(end - start,2)
    print(f'\nSniffles finished compressing and indexing in {runtime} seconds')
    start = time.time()

    #command list for generating consensus fasta
    cmds = []
    for id in ids:
        #use bcftools to get consensus fasta
        cmd = f'{libPath}/bin/bcftools consensus -f {reference} {outDir}/consensus/{id}.vcf.gz -o {outDir}/consensus/{id}.fasta'
        cmds.append(cmd)
        #add id as finished on tracking
        readData.addData('consensus',id,f'{outDir}/consensus/{id}.fasta')
    #start multiprocessing
    pool.starmap(proc, [[runCFG,i] for i in cmds])

    #determine runtime of processes
    end = time.time()
    runtime = round(end - overall_start,2)
    print(f'\nSniffles finished generating consensus sequence in {runtime} seconds')

    #if mapping reads to consensus is specified run mapping
    if runCFG['exec']['mapToConsensus']:
        refs = []
        for id in ids:
            refs.append([id,readData.data['consensus'][id]])
        mapping(readData,runCFG,threads,ids,refs=refs,jobtype='map-consensus')

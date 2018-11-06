import subprocess as sub
import os
import shlex
from mapping import mapping
import multiprocessing as mp
import time
from sniffProc import proc,init
from sc import procTitle,checkexists

def removeDuplicates(readData,runCFG,threads='1',ids=''):
    #inital parameters
    libPath = runCFG['libPath']
    outDir =runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])

    #if no id supplied get list of all ids
    if not ids:
        ids = readData.idList

    #notify starting to remove duplicates
    procTitle('Remove Duplicates')
    print('\nSniffles: Removing duplicate reads')
    #get time at start
    start = time.time()

    #generate commands
    cmds = []
    sort_cmds = []
    checkexists(os.path.join(outDir,'nodups'))
    for id in ids:
        #determine which samfile to use if it has been normalized
        if id in readData.data['normalized']:
            samfile = f'{outDir}/normalized_mapping/{id}.sam'
            cmd = f'{libPath}/bin/samtools view -b {samfile} | samtools sort | samtools view -h > {outDir}/nodups/{id}_pre.sam'
        else:
            samfile = f'{outDir}/inital_mapping/{id}.sam'
            cmd = f'{libPath}/bin/samtools view -b {samfile} | samtools sort | samtools view -h > {outDir}/nodups/{id}_pre.sam'
        sort_cmds.append(cmd)

        #remove duplicate reads command
        cmd = f'java -Xmx2g -jar {libPath}/picard/picard.jar MarkDuplicates I={outDir}/nodups/{id}_pre.sam O={outDir}/nodups/{id}.sam REMOVE_DUPLICATES=true M={id}.removeDupMetrics.txt'
        cmds.append(cmd)

    #set up multiprocessing
    #start multiprocessing
    lock = mp.Lock()
    pool = mp.Pool(processes=threads,initializer=init,initargs=(lock,))

    #denote start of remove duplicate reads in logs
    with open(logfile,'a') as outlog:
        outlog.write('**********************\n')
        outlog.write('Remove Duplicate Reads\n')
    #begin multiprocessing
    pool.starmap(proc,[[runCFG,i,'','',True] for i in sort_cmds])
    pool.starmap(proc, [[runCFG,i] for i in cmds])
    #get time at end
    end = time.time()
    with open(logfile,'a') as outlog:
        outlog.write('**********************\n')
    #determine runtime of processes
    runtime = round(end - start,2)
    print(f'\nSniffles finished removing duplicates in {runtime} seconds')

    #cleanup
    for id in ids:
        #remove intermediate files
        os.remove(outDir+'/nodups/'+'{id}_pre.sam'.format(id=id))
        #add id to finished list
        readData.addData('rmDuplicates',id,f'{outDir}/nodups/{id}.sam')

def normCoverage(readData,runCFG,threads='1',ids=''):
    #inital parameters
    libPath = runCFG['libPath']
    outDir =runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])

    #if no id supplied get list of all ids
    if not ids:
        ids = readData.idList

    #generate commands
    format_cmds = []
    norm_cmds = []
    checkexists(os.path.join(outDir,'normalized'))
    for id in ids:
        #get reads from mapped samfile
        cmd = f'{libPath}/bin/samtools fastq {outDir}/inital_mapping/{id}.sam -1 {outDir}/normalized/{id}_1.fastq -2 {outDir}/normalized/{id}_2.fastq'
        format_cmds.append(cmd)

        #run bbnorm
        cov = runCFG['exec']['coverageNormDepth']
        cmd = f'{libPath}/bbmap/bbnorm.sh in={outDir}/normalized/{id}_1.fastq in2={outDir}/normalized/{id}_2.fastq out={outDir}/normalized/{id}_normalized.fastq target={cov} threads={threads}'
        norm_cmds.append(cmd)

    #notify starting to remove duplicates
    procTitle('Normalize Coverage')
    print('\nSniffles: Normalizing read coverage')

    #denote start of remove duplicate reads in logs
    with open(logfile,'a') as outlog:
        outlog.write('********************\n')
        outlog.write('Normalizing coverage\n')

    #get time at start
    start = time.time()

    #set up multiprocessing
    #start multiprocessing
    lock = mp.Lock()

    pool = mp.Pool(processes=threads,initializer=init,initargs=(lock,))
    pool.starmap(proc, [[runCFG,i] for i in format_cmds])

    #NOTE: although normalizing is setup for multiprocessing bbnorm
    #uses all available memory, thus can only be run serialy
    pool = mp.Pool(processes=1,initializer=init,initargs=(lock,))
    pool.starmap(proc, [[runCFG,i] for i in norm_cmds])

    #get time at end
    end = time.time()
    with open(logfile,'a') as outlog:
        outlog.write('********************\n')

    #determine runtime of processes
    runtime = round(end - start,2)
    print(f'\nSniffles finished normalizing read coverage in {runtime} seconds')

    #cleanup
    for id in ids:
        #os.remove(f'{outDir}/normalized/{id}_1.fastq')
        #os.remove(f'{outDir}/normalized/{id}_2.fastq')
        #add to tracker
        readData.addData('normalized',id,f'{outDir}/normalized/{id}_normalized.fastq')

    mapping(readData,runCFG,threads,ids,jobtype='map-normalized')

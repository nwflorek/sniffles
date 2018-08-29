import subprocess as sub
import os
import shlex
from mapping import mapping
import multiprocessing as mp
import time
from sniffProc import proc,init

def removeDuplicates(readData,runCFG,threads='1',ids=''):
    #inital parameters
    libPath = runCFG['libPath']
    outDir =runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])

    #if no id supplied get list of all ids
    if not ids:
        ids = readData.idList

    #notify starting to remove duplicates
    print('\n********************************')
    print('Sniffles: Removing duplicate reads')
    #get time at start
    start = time.time()

    #generate commands
    cmds = []
    for id in ids:
        #sort samfile
        cmd01 = f'{libPath}/bin/samtools view -b {id}.sam'
        cmd01 = shlex.split(cmd01)
        cmd02 = f'{libPath}/bin/samtools sort'
        cmd02 = shlex.split(cmd02)
        cmd03 = f'{libPath}/bin/samtools view -h'
        cmd03 = shlex.split(cmd03)
        with open(outDir+'/'+'{id}_sorted.sam'.format(id=id),'w') as outsam:
            c1 = sub.Popen(cmd01,stdout=sub.PIPE,cwd=outDir)
            c2 = sub.Popen(cmd02,stdin=c1.stdout,stdout=sub.PIPE,cwd=outDir)
            c3 = sub.Popen(cmd03,stdin=c2.stdout,stdout=outsam,cwd=outDir)
            c3.wait()

        #remove duplicate reads command
        cmd = f'java -jar {libPath}/picard/picard.jar MarkDuplicates I={id}_sorted.sam O={id}_nodups.sam REMOVE_DUPLICATES=true M={id}.removeDupMetrics.txt'
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
    pool.starmap(proc, [[runCFG,i] for i in cmds])
    #get time at end
    end = time.time()
    with open(logfile,'a') as outlog:
        outlog.write('**********************\n')
    #determine runtime of processes
    runtime = end - start
    print(f'\nSniffles finished removing duplicates in {runtime} seconds')
    print('\n********************************')

    #cleanup
    for id in ids:
        #remove intermediate files
        os.remove(outDir+'/'+'{id}_sorted.sam'.format(id=id))
        os.remove(outDir+'/'+'{id}.sam'.format(id=id))
        #add id to finished list
        readData.addData('rmDuplicates',id,f'{id}_nodups.sam')

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
    for id in ids:
        #determine which samfile to use if duplicates have been removed
        if id in readData.data['rmDuplicates']:
            samfile = f'{id}_nodups.sam'
        else:
            samfile = f'{id}.sam'

        #get reads from samfile
        cmd = f'{libPath}/bbmap/reformat.sh in={samfile} out={id}_adjusted.fastq'
        format_cmds.append(cmd)

        #run bbnorm
        cov = runCFG['exec']['coverageNormDepth']
        cmd = f'{libPath}/bbmap/bbnorm.sh in={id}_adjusted.fastq out={id}_normalized.fastq target={cov}'
        norm_cmds.append(cmd)

    #set up multiprocessing
    #start multiprocessing
    lock = mp.Lock()
    
    #NOTE: although normalizing is setup for multiprocessing bbnorm
    #uses all available memory, thus can only be run serialy
    pool = mp.Pool(processes=1,initializer=init,initargs=(lock,))
    #notify starting to remove duplicates
    print('\n***********************************')
    print('\nSniffles: Normalizing read coverage')

    #denote start of remove duplicate reads in logs
    with open(logfile,'a') as outlog:
        outlog.write('********************\n')
        outlog.write('Normalizing coverage\n')

    #get time at start
    start = time.time()

    #begin multiprocessing
    pool.starmap(proc, [[runCFG,i] for i in format_cmds])
    pool.starmap(proc, [[runCFG,i] for i in norm_cmds])

    #get time at end
    end = time.time()
    with open(logfile,'a') as outlog:
        outlog.write('********************\n')

    #determine runtime of processes
    runtime = end - start
    print(f'\nSniffles finished normalizing read coverage in {runtime} seconds')
    print('\n***********************************')

    #cleanup
    for id in ids:
        os.remove(f'{outDir}/{id}_adjusted.fastq')
        #add to tracker
        readData.addData('normalized',id,f'{outDir}/{id}_normalized.fastq')

    mapping(readData,runCFG,threads,ids)

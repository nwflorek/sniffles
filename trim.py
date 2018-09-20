import os
import shlex
import subprocess as sub
import multiprocessing as mp
import time
from sniffProc import proc,init
from sc import procTitle

def trimmomatic(readData,runCFG,threads,ids=''):

    #parameters
    minlength = runCFG['trimmomatic']['minlength']
    windowsize = runCFG['trimmomatic']['windowSize']
    qscore = runCFG['trimmomatic']['qscore']
    adapterpath = runCFG['trimmomatic']['adaptersPath']
    outDir = runCFG['exec']['outdir']
    libPath=runCFG['libPath']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])

    #set up list of ids to trim
    if not ids:
        ids = readData.idList

    #generate commands for each trim job
    cmds = []
    for id in ids:
        #main command
        main_cmd = f'java -jar {libPath}/trimmomatic/trimmomatic-0.36.jar '

        #get read path
        if readData.rawdata[id]:
            read1 = readData.rawdata[id][0]
            read2 = readData.rawdata[id][1]
        #determine args
        if runCFG['trimmomatic']['removeAdapters']:
            if runCFG['trimmomatic']['paired']:
                args = f'PE {read1} {read2} -baseout {id}_trimmed.fastq.gz ILLUMINACLIP:{adapterpath}:1:30:10 SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
                readData.addData('trimmed',id,f'{outDir}/trimmed/{id}_trimmed_1P.fastq.gz',f'trimmed/{id}_trimmed_2P.fastq.gz')
            else:
                args = f'SE {read1} -baseout {id}_trimmed.fastq.gz ILLUMINACLIP:{adapterpath}:1:30:10 SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
                readData.addData('trimmed',id,f'{outDir}/trimmed/{id}_trimmed.fastq.gz')
        else:
            if runCFG['trimmomatic']['paired']:
                args = f'PE {read1} {read2} -baseout {id}_trimmed.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
                readData.addData('trimmed',id,f'{outDir}/trimmed/{id}_trimmed_1P.fastq.gz',f'trimmed/{id}_trimmed_2P.fastq.gz')
            else:
                args = f'SE {read1} -baseout {id}_trimmed.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
                readData.addData('trimmed',id,f'{outDir}/trimmed/{id}_trimmed.fastq.gz')

        #prepare command and add to list
        sample_cmd = main_cmd+args
        cmds.append(sample_cmd)

    #make out dir if it doesn't already exist
    try:
        os.mkdir(os.path.join(outDir,'trimmed'))
    except:
        pass

    #set up multiprocessing
    #start multiprocessing
    lock = mp.Lock()
    pool = mp.Pool(processes=threads,initializer=init,initargs=(lock,))
    #notify starting trimming
    procTitle('Quality Trimming')
    print('\nSniffles: Started quality trimming')
    #start timer
    start = time.time()
    #denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Trimmomatic\n')
    #begin multiprocessing
    pool.starmap(proc, [[runCFG,i,'trimmed'] for i in cmds])
    #denote end of logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
    #get time
    end = time.time()
    #determine runtime of processes
    runtime = round(end - start,2)
    print(f'\nFinished trimming in {runtime} seconds')

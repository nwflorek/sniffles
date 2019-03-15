import os
import shlex
import subprocess as sub
import multiprocessing as mp
import time
from sc import procTitle
import calldocker as cd

def trim(readData,runCFG,threads,ids=''):

    #parameters
    minlength = runCFG['trimmomatic']['minlength']
    windowsize = runCFG['trimmomatic']['windowSize']
    qscore = runCFG['trimmomatic']['qscore']
    adapterpath = "/tools/adapters/" + runCFG['trimmomatic']['adaptersFileName']
    outDir = runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])

    #set up list of ids to trim
    if not ids:
        ids = readData.ids

    #generate commands for each trim job
    cmds = []
    for id in ids:
        #main command
        main_cmd = f'java -jar /tools/trimmomatic.jar '

        #get read path
        if readData.reads[id]:
            read_path = os.path.dirname(os.path.abspath(readData.reads[id].fwd))
            read1_basename = os.path.basename(readData.reads[id].fwd)
            read2_basename = os.path.basename(readData.reads[id].rev)

        #determine args
        if runCFG['trimmomatic']['removeAdapters']:
            if runCFG['trimmomatic']['paired']:
                args = f'PE {read1_basename} {read2_basename} -baseout /output/{id}_trimmed.fastq.gz ILLUMINACLIP:{adapterpath}:1:30:10 SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
                readData.add_runtime('trimmed',id,f'{outDir}/trimmed/{id}_trimmed_1P.fastq.gz',f'{outDir}/trimmed/{id}_trimmed_2P.fastq.gz')
            else:
                args = f'SE {read1_basename} -baseout /output/{id}_trimmed.fastq.gz ILLUMINACLIP:{adapterpath}:1:30:10 SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
                readData.add_runtime('trimmed',id,f'{outDir}/trimmed/{id}_trimmed.fastq.gz')
        else:
            if runCFG['trimmomatic']['paired']:
                args = f'PE {read1_basename} {read2_basename} -baseout /output/{id}_trimmed.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
                readData.add_runtime('trimmed',id,f'{outDir}/trimmed/{id}_trimmed_1P.fastq.gz',f'{outDir}/trimmed/{id}_trimmed_2P.fastq.gz')
            else:
                args = f'SE {read1_basename} -baseout /output/{id}_trimmed.fastq.gz SLIDINGWINDOW:{windowsize}:{qscore} MINLEN:{minlength}'
                readData.add_runtime('trimmed',id,f'{outDir}/trimmed/{id}_trimmed.fastq.gz')

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
    pool = mp.Pool(processes=threads)
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
        results = pool.starmap_async(cd.call,[[cmd,'/data',{read_path:"/data",os.path.join(outDir,'trimmed'):"/output"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')
    #get time
    end = time.time()
    #determine runtime of processes
    runtime = round(end - start,2)
    print(f'\nFinished trimming in {runtime} seconds')

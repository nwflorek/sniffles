import os
import shlex
import subprocess as sub
import multiprocessing as mp

#define process for multiprocessing
def trim(outDir,cmd):
    cmd = shlex.split(cmd)
    p = sub.Popen(cmd,cwd=os.path.join(outDir,'trimmed'),stdout=sub.PIPE)
    lock.acquire()
    for line in p.stdout:
        print(line)
    lock.release()
    p.wait()
#make a lock global for child workers
def init(l):
    global lock
    lock = l

def trimmomatic(readData,runCFG,threads,ids=''):

    #parameters
    minlength = runCFG['trimmomatic']['minlength']
    windowsize = runCFG['trimmomatic']['windowSize']
    qscore = runCFG['trimmomatic']['qscore']
    adapterpath = runCFG['trimmomatic']['adaptersPath']
    outDir = runCFG['exec']['outdir']
    libPath=runCFG['libPath']

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
    pool.starmap(trim, [[outDir,i] for i in cmds])

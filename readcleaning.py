import subprocess as sub
import os
import shlex
import calldocker as cd
import multiprocessing as mp
import time
from sc import procTitle,checkexists

def removeDuplicates(runCFG,bam_files,threads='1'):
    #inital parameters
    outDir =runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])
    checkexists(os.path.join(outDir,'rm_dups'))
    outDir = os.path.join(outDir,'rm_dups')

    #notify starting to remove duplicates
    procTitle('Remove Duplicates')
    print('\nSniffles: Removing duplicate reads')
    #get time at start
    start = time.time()

    #generate commands
    cmds = []
    output_list = []
    for path in bam_files:
        full_path = os.path.abspath(path)
        file_name = os.path.basename(full_path)
        path = os.path.dirname(full_path)
        id = file_name.split(".")[0]

        #remove duplicate reads command
        cmd = f'java -Xmx2g -jar /tools/picard.jar MarkDuplicates I=/in_dir/{id}.bam O=/out_dir/{id}.bam REMOVE_DUPLICATES=true M=/out_dir/{id}.removeDupMetrics.txt'
        cmds.append(cmd)

        #add id to finished list
        output_list.append(os.path.join(outDir,f'{id}.bam'))

    #set up multiprocessing
    #start multiprocessing
    pool = mp.Pool(processes=threads)

    #denote start of remove duplicate reads in logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Removing Duplicates\n')
        #start multiprocessing
        results = pool.starmap_async(cd.call,[[cmd,'/reads',{path:"/in_dir",outDir:"/out_dir"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')
    #get time at end
    end = time.time()

    #determine runtime of processes
    runtime = round(end - start,2)
    print(f'\nSniffles: Finished removing duplicates in {runtime} seconds')
    return output_list


def normCoverage(runCFG,bam_files,threads='1'):
    #NOTE: normalizing with bbnorm uses all available memory, thus can only be run serialy

    #inital parameters
    outDir =runCFG['exec']['outdir']
    checkexists(os.path.join(outDir,'normalized'))
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])
    outDir = os.path.join(outDir,'normalized')

    #notify starting to remove duplicates
    procTitle('Normalize Coverage')
    print('\nSniffles: Normalizing read coverage')

    #get time at start
    start = time.time()

    #denote start of remove duplicate reads in logs
    with open(logfile,'a') as outlog:
        outlog.write('********************\n')
        outlog.write('Normalizing coverage\n')

        #run normalization
        output_list = []
        for path in bam_files:
            full_path = os.path.abspath(path)
            file_name = os.path.basename(full_path)
            path = os.path.dirname(full_path)
            id = file_name.split('.')[0]

            #get reads from mapped bamfile
            cmd_get_reads = f'bash -c \'samtools fastq /bam_files/{id}.bam -1 /out_dir/{id}_mapped_1.fastq -2 /out_dir/{id}_mapped_2.fastq && '

            #run seqtk to subsample reads
            total_reads = runCFG['exec']['totalReads']
            cmd_normalization = f'seqtk sample -s100 /out_dir/{id}_mapped_1.fastq {total_reads} > {id}_1.fastq && seqtk sample -s100 /out_dir/{id}_mapped_2.fastq {total_reads} > {id}_2.fastq\''

            #start docker containers and run
            outlog.write(f'{id}-----------\n')
            stdout=cd.call(cmd_get_reads+cmd_normalization,'/out_dir',{path:"/bam_files",outDir:"/out_dir"})
            outlog.write(stdout)
            outlog.write(f'-----------\n')

            output_list.append([os.path.join(outDir,f'{id}_1.fastq'),os.path.join(outDir,f'{id}_2.fastq')])


            #cleanup
            try:
                os.remove(f'{outDir}/{id}_mapped_1.fastq')
            except:
                pass
            try:
                os.remove(f'{outDir}/{id}_mapped_2.fastq')
            except:
                pass
        outlog.write('********************\n')

    #get time at end
    end = time.time()

    #determine runtime of processes
    runtime = round(end - start,2)
    print(f'\nSniffles: Finished normalizing read coverage in {runtime} seconds')
    return output_list

import os
import shlex
import subprocess as sub
import multiprocessing as mp
import time
from sniffProc import proc,init
from sc import procTitle

def mapping(readData,runCFG,threads='1',ids='',refs=None):
    #inital parameters
    if not refs:
        reference_sequence = os.path.abspath(runCFG['exec']['referenceSequence'])
        reference_sequence_name = os.path.basename(reference_sequence)
    libPath = runCFG['libPath']
    outDir = runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])

    #use id if provided otherwise get list
    if not ids:
        ids = readData.idList

    #index reference
    if 'bowtieindexed' not in readData.data:
        cmd = f'{libPath}/bin/bowtie2-build {reference_sequence} {reference_sequence_name}'
        cmd = shlex.split(cmd)
        with open(logfile,'a') as outlog:
            outlog.write("*************************\n")
            outlog.write("Bowtie indexing reference\n")
        with open(logfile,'a') as outlog:
            sub.Popen(cmd,cwd=outDir,stdout=outlog,stderr=outlog).wait()
            outlog.write("*************************\n")
        readData.data['bowtieindexed'] = True
    ref_dict = {}
    if refs:
        for ref in refs:
            reference_sequence = os.path.abspath(ref[1][0])
            reference_sequence_name = os.path.basename(ref[1][0])
            cmd = f'{libPath}/bin/bowtie2-build {reference_sequence} {reference_sequence_name}'
            cmd = shlex.split(cmd)
            with open(logfile,'a') as outlog:
                outlog.write("*************************\n")
                outlog.write("Bowtie indexing reference\n")
            with open(logfile,'a') as outlog:
                sub.Popen(cmd,cwd=outDir,stdout=outlog,stderr=outlog).wait()
                outlog.write("*************************\n")
            readData.data['bowtieindexed'] = True
            ref_dict[ref[0]] = reference_sequence_name

    #generate mapping commands
    cmds = []
    for id in ids:
        #determine reads
        if 'trimmed' in readData.data and id not in readData.data['mapProgress']['trimmed']:
            #TODO add status for unpaired read information
            reads = [readData.data['trimmed'][id][0],readData.data['trimmed'][id][1]]
            readData.data['mapProgress']['trimmed'].append(id)
            interleaved = False

        elif 'normalized' in readData.data and id not in readData.data['mapProgress']['normalized']:
            reads = readData.data['normalized'][id]
            readData.data['mapProgress']['normalized'].append(id)
            interleaved = True

        elif 'consensus' in readData.data and id not in readData.data['mapProgress']['consensus']:
            reads = [readData.data['trimmed'][id][0],readData.data['trimmed'][id][1]]
            readData.data['mapProgress']['normalized'].append(id)
            interleaved = False
        else:
            print(f'There was an error mapping {id}.')
            exit()

        #determine interleaved or not
        #interleaved cmd
        if ref_dict:
            reference_sequence_name = ref_dict[id]
            sam = 'consensus'
        else:
            sam = 'remapped'
        if interleaved:
            interleaved_cmd = f"{libPath}/bin/bowtie2 -x {reference_sequence_name} --interleaved {reads[0]} -S {id}_{sam}.sam -p 2 --local"
            cmds.append(interleaved_cmd)
        #split cmd
        else:
            split_cmd = f"{libPath}/bin/bowtie2 -x {reference_sequence_name} -1 {reads[0]} -2 {reads[1]} -S {id}.sam -p 2 --local"
            cmds.append(split_cmd)

    #set up multiprocessing
    #start multiprocessing
    lock = mp.Lock()
    pool = mp.Pool(processes=1,initializer=init,initargs=(lock,))
    #notify starting mapping
    procTitle('Read Mapping')
    print('\nSniffles: Started mapping')
    #denote start of mapping in logs
    with open(logfile,'a') as outlog:
        outlog.write('*******\n')
        outlog.write('Mapping\n')
    #get start time
    start = time.time()
    #start multiprocessing
    pool.starmap(proc, [[runCFG,i] for i in cmds])
    #get end time
    end = time.time()
    #denote end of mapping in log
    with open(logfile,'a') as outlog:
        outlog.write('*******\n')
    #get total runtime
    runtime = end - start
    print(f'\nSniffles finished mapping in {runtime} seconds')

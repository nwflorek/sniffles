import os
import shlex
import subprocess as sub
import multiprocessing as mp
import time
from sniffProc import proc,init
from sc import procTitle,checkexists
from readdepth import average_depth

def mapping(readData,runCFG,threads='1',ids='',refs=None,jobtype=None):
    #TODO update this to a better method
    #determine how many concurrent processes to run at a time
    concurrency = int(threads / 2)

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
        checkexists(f'{outDir}/indexes')
        for ref in refs:
            ref_data = ref[1][0]
            reference_sequence = ref_data
            reference_sequence_name = os.path.basename(ref_data)
            cmd = f'{libPath}/bin/bowtie2-build {reference_sequence} {reference_sequence_name}'
            cmd = shlex.split(cmd)
            with open(logfile,'a') as outlog:
                outlog.write("*************************\n")
                outlog.write("Bowtie indexing reference\n")
            with open(logfile,'a') as outlog:
                sub.Popen(cmd,cwd=os.path.join(outDir,'indexes'),stdout=outlog,stderr=outlog).wait()
                outlog.write("*************************\n")
            readData.data['bowtieindexed'] = True
            ref_dict[ref[0]] = os.path.join(outDir,'indexes',reference_sequence_name)

    #generate mapping commands
    cmds = []
    check_depth = False
    for id in ids:
        #setup for multiple references
        if ref_dict:
            reference_sequence_name = ref_dict[id]

        #determine jobtype and generate mapping command
        if 'map-trimmed' == jobtype:
            #TODO add status for unpaired read information
            #get read paths from readdata object (stored by trim)
            reads = [readData.data['trimmed'][id][0],readData.data['trimmed'][id][1]]
            #add id to progress tracker
            readData.data['mapProgress']['trimmed'].append(id)
            #check output folder exists
            checkexists(os.path.join(outDir,'inital_mapping'))
            #generate command
            cmd = f"{libPath}/bin/bowtie2 -x {reference_sequence_name} -1 {reads[0]} -2 {reads[1]} -p 2 --local | {libPath}/bin/samtools view -bS - | {libPath}/bin/samtools sort > {outDir}/inital_mapping/{id}.bam"
            cmds.append(cmd)
            check_depth = True

        elif 'map-normalized' == jobtype:
            #get read path from tracker
            reads = readData.data['normalized'][id]
            #add ids to mapping tracker
            readData.data['mapProgress']['normalized'].append(id)
            #check output folder exists
            checkexists(os.path.join(outDir,'normalized_mapping'))
            cmd = f"{libPath}/bin/bowtie2 -x {reference_sequence_name} --interleaved {reads[0]} -p 2 --local | {libPath}/bin/samtools view -bS - | {libPath}/bin/samtools sort > {outDir}/normalized_mapping/{id}.bam"
            cmds.append(cmd)

        elif 'map-consensus' == jobtype:
            #get reads path from tracker
            reads = [readData.data['trimmed'][id][0],readData.data['trimmed'][id][1]]
            #add ids to mapping tracker
            readData.data['mapProgress']['normalized'].append(id)
            #check output folder exists
            checkexists(os.path.join(outDir,'consensus_mapping'))
            #generate command
            cmd = f"{libPath}/bin/bowtie2 -x {reference_sequence_name} -1 {reads[0]} -2 {reads[1]} -p 2 --local | {libPath}/bin/samtools view -bS - | {libPath}/bin/samtools sort > {outDir}/consensus_mapping/{id}.bam"
            cmds.append(cmd)
        else:
            print(f'There was an error mapping {id}.')
            exit()

    #set up multiprocessing
    #start multiprocessing
    lock = mp.Lock()
    pool = mp.Pool(processes=concurrency,initializer=init,initargs=(lock,))
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
    pool.starmap(proc, [[runCFG,i,'','',True] for i in cmds])
    #determine average depth for each isolate
    if check_depth:
        passingDepth = []
        for id in ids:
            depth = average_depth(f'{outDir}/inital_mapping/{id}.bam')
            readData.data['mapData']['avgDepth'][id] = depth
            #only keep isolates that pass the minimumaveragedepth
            if int(float(depth)) >= runCFG['exec']['minimumAverageDepth']:
                passingDepth.append(id)
        readData.idList = passingDepth

    #get end time
    end = time.time()
    #denote end of mapping in log
    with open(logfile,'a') as outlog:
        outlog.write('*******\n')
    #get total runtime
    runtime = round(end - start,2)
    print(f'\nSniffles finished mapping in {runtime} seconds')

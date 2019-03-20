import os
import shlex
import subprocess as sub
import multiprocessing as mp
import time
from sc import procTitle,checkexists,cpu_count
import calldocker as cd
from readdepth import average_depth
from shutil import copyfile

def indexing(runCFG,*paths):
    logfile = os.path.join(runCFG['exec']['outdir'],runCFG['exec']['logfile'])
    outDir = runCFG['exec']['outdir'] + '/ref_sequence'
    os.mkdir(outDir)
    procTitle('Indexing Reference Genome')
    for path in paths:
        reference_sequence_abspath = os.path.abspath(path)
        reference_sequence_name = os.path.basename(reference_sequence_abspath)
        #index reference
        cmd = f'bowtie2-build {reference_sequence_name} {reference_sequence_name}'
        with open(logfile,'a') as outlog:
            outlog.write("*************************\n")
            outlog.write("Bowtie2 indexing the reference\n")
            copyfile(reference_sequence_abspath,os.path.join(outDir,reference_sequence_name))
            outlog.write(cd.call(cmd,'/data',{outDir:"/data"}))
            outlog.write("*************************\n")

#mapping parameters
#param_path - list of paths in the following order
#   [(id,path to fwd read, path to rev read, path to reference),...]
#   if not paired end then "path to rev read" will be empty
#threads - number of threads for bowtie2
#outDir - the output directory
def mapping(runCFG,param_paths,outDir,threads='1'):

    logfile = os.path.join(runCFG['exec']['outdir'],runCFG['exec']['logfile'])

    num_jobs,num_threads = cpu_count(threads)

    cmds = []
    read_path = ''
    ref_path = ''
    for param_path in param_paths:
        id = param_path[0]
        read1 = os.path.basename(param_path[1])
        read2 = os.path.basename(param_path[2])
        read_path = os.path.dirname(os.path.abspath(param_path[1]))
        ref_path = runCFG['exec']['outdir'] + '/ref_sequence'
        reference_sequence_name = os.path.basename(param_path[3])

        #check output folder exists
        checkexists(os.path.join(outDir))
        #generate command
        cmd = f"bash -c \'bowtie2 -x {reference_sequence_name} -1 /reads/{read1} -2 /reads/{read2} -p {num_threads} --local | samtools view -bS | samtools sort -o /output/{id}.bam\'"
        cmds.append(cmd)

    #set up multiprocessing
    #start multiprocessing
    pool = mp.Pool(processes=num_jobs)
    #notify starting mapping
    procTitle('Mapping Reads')
    print('\nSniffles: Started mapping')
    #get start time
    start = time.time()
    #denote start of mapping in logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Mapping\n')
        #start multiprocessing
        results = pool.starmap_async(cd.call,[[cmd,'/reads',{ref_path:"/reference",read_path:"/reads",outDir:"/output"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')
    #get end time
    end = time.time()
    #get total runtime
    runtime = round(end - start,2)
    print(f'\nSniffles finished mapping in {runtime} seconds')
'''
    #determine average depth for each isolate
    if check_depth:
        passingDepth = []
        for id in ids:
            depth = average_depth(f'{outDir}/inital_mapping/{id}.bam')
            if 'error' in depth:
                readData.data['mapData']['avgDepth'][id] = depth
                continue
            readData.data['mapData']['avgDepth'][id] = depth
            #only keep isolates that pass the minimumaveragedepth
            if int(float(depth)) >= runCFG['exec']['minimumAverageDepth']:
                passingDepth.append(id)
        readData.idList = passingDepth
'''

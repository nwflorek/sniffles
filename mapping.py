import os, re
import shlex
import subprocess as sub
import multiprocessing as mp
import time
from sc import procTitle,checkexists,cpu_count
import calldocker as cd
from shutil import copyfile

def indexing(runCFG,*paths):
    logfile = os.path.join(runCFG['exec']['outdir'],runCFG['exec']['logfile'])
    outDir = runCFG['exec']['outdir'] + '/ref_sequence'
    checkexists(outDir)
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
    output_bam_list = []
    for param_path in param_paths:
        id = param_path[0]
        read1 = os.path.basename(param_path[1])
        read2 = os.path.basename(param_path[2])
        read_path = os.path.dirname(os.path.abspath(param_path[1]))
        ref_path = runCFG['exec']['outdir'] + '/ref_sequence'
        reference_sequence_name = os.path.basename(param_path[3])

        #check output folder exists
        checkexists(os.path.join(outDir))
        if read2 != '':
            #generate command for paired end
            cmd = f"bash -c \'bowtie2 -x {reference_sequence_name} -1 /reads/{read1} -2 /reads/{read2} -p {num_threads} --local | samtools view -bS | samtools sort -o /output/{id}.bam\'"
        else:
            #generate command for interleaved
            cmd = f"bash -c \'bowtie2 -x {reference_sequence_name} --interleaved /reads/{read1} -p {num_threads} --local | samtools view -bS | samtools sort -o /output/{id}.bam\'"
        cmds.append(cmd)

        #data for next stage
        output_bam_list.append(os.path.join(outDir,f'{id}.bam'))

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
    print(f'\nSniffles: Finished mapping in {runtime} seconds')
    return output_bam_list

def average_depth(runCFG,bam_list,inDir,outDir):
    #check that output folder exists
    checkexists(os.path.join(outDir))
    #setup inital parameters
    ref_path = runCFG['exec']['outdir'] + '/ref_sequence'
    reference_sequence_name = os.path.basename(runCFG['exec']['referenceSequence'])
    logfile = os.path.join(runCFG['exec']['outdir'],runCFG['exec']['logfile'])

    #bam file list that will meet threashold
    filtered_bam_list = []

    #create logfile that will hold all average depths
    with open(os.path.join(outDir,'average_depth.log'),'w') as outdepth:

        #loop through each bam
        for bam in bam_list:
            filename = os.path.basename(bam)
            id = filename.split('.')[0]
            #open the log to log output
            with open(logfile,'a') as outlog:
                outlog.write('***********\n')
                outlog.write('Coverage\n')
                #generate command for the current bam file
                cmd = f'bash -c "samtools view /indata/{filename} > {id}.tmp.sam  && /tools/bbmap/pileup.sh in={id}.tmp.sam out={id}_coverage.csv ref=/reference/{reference_sequence_name} && rm {id}.tmp.sam"'
                #use docker to run the command
                output = cd.call(cmd,'/outdata',{inDir:"/indata",outDir:"/outdata",ref_path:"/reference"})
                #record the output in the log
                outlog.write(output)
                #denote end of logs
                outlog.write('***********\n')

            #only add isolates that pass average depth and percent of reference covered
            percent_cov = 0
            avg_cov = 0
            #parse lines of stdout for info we need
            for line in output.splitlines():
                if "Percent of reference bases covered:" in line:
                    match = re.search('[0-9,.]+',line)
                    if match:
                        percent_cov = float(match[0])

                if "Average coverage:" in line:
                    match = re.search('[0-9,.]+',line)
                    if match:
                        avg_cov = float(match[0])
            #check against the config for min thresholds
            if avg_cov >= runCFG['exec']['minimumAverageDepth'] and percent_cov >= runCFG['exec']['percentRefCovered']:
                #record
                outdepth.write(f'{id},{avg_cov},{percent_cov},Pass\n')
                filtered_bam_list.append(bam)
            else:
                outdepth.write(f'{id},{avg_cov},{percent_cov},Fail\n')
    #return only bam files that meet coverage requirements
    return filtered_bam_list

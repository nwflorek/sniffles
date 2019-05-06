import os
import shlex
import subprocess as sub
import time
from sc import procTitle,checkexists
import multiprocessing as mp
import calldocker as cd

def consensus(runCFG,bam_list,threads='1'):
    #inital parameters
    outDir =runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])
    outDir = os.path.join(outDir,'consensus')
    checkexists(outDir)

    #notify starting mapping
    procTitle('Generate Consensus')
    print('\nSniffles: Started generating consensus vcf')

    #get start time
    overall_start = time.time()
    start = time.time()

    #set reference sequence
    reference_sequence_abspath = os.path.abspath(runCFG['exec']['referenceSequence'])
    reference_sequence_name = os.path.basename(reference_sequence_abspath)
    reference_sequence_dir = runCFG['exec']['outdir'] + '/ref_sequence'

    #command list
    cmds = []
    vcf_list = []
    for path in bam_list:
        full_path = os.path.abspath(path)
        file_name = os.path.basename(full_path)
        path = os.path.dirname(full_path)
        id = file_name.split(".")[0]
        #run varscan mpileup2cns to generate vcf with consensus information
        minCov = runCFG['snpcalling']['minCoverage']
        quality = runCFG['snpcalling']['snpQualityThreshold']
        freq = runCFG['snpcalling']['consensusFrequency']

        #make multiway pileup using samtools
        cmd1 = f'bash -c \'samtools mpileup -d 1000000 /infile/{file_name} -f /ref/{reference_sequence_name} -o {id}.pileup && '

        cmd2 = f'java -jar /tools/varscan.jar mpileup2cns {id}.pileup --min-coverage {minCov} --min-avg-qual {quality} --min-var-freq {freq} --strand-filter 1 --output-vcf 1 > {id}.vcf\''
        cmds.append(cmd1 + cmd2)
        vcf_list.append(os.path.join(outDir,f'{id}.vcf'))

    #setup multiprocessing
    pool = mp.Pool(processes=threads)

    #start multiprocessing
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Generating Consensus\n')
        #start multiprocessing
        results = pool.starmap_async(cd.call,[[cmd,'/outfile',{reference_sequence_dir:"/ref",path:"/infile",outDir:"/outfile"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')

    #check if vcf file is empty, if it is skip id and remove vcf file
    filtered_vcf_list = []
    for path in vcf_list:
        try:
            if os.path.getsize(path)>0:
                filtered_vcf_list.append(path)
            else:
                os.remove(path)
        except:
            pass

    end = time.time()
    runtime = round(end - start,2)
    print(f'\nSniffles: Finished generating the consensus vcf in {runtime} seconds')
    start = time.time()
    print(f'\nSniffles: Generating consensus fasta')

    #command list for compressing files
    cmds = []
    out_fasta = []
    for vcf in filtered_vcf_list:
        full_path = os.path.abspath(vcf)
        file_name = os.path.basename(full_path)
        path = os.path.dirname(full_path)
        id = file_name.split(".")[0]
        #compress vcf file with bgzip
        cmd = f'bash -c \'bgzip {id}.vcf && tabix {id}.vcf.gz && bcftools consensus -f /ref/{reference_sequence_name} {id}.vcf.gz -o {id}.fasta\''
        out_fasta.append(os.path.join(outDir,f'{id}.fasta'))
        cmds.append(cmd)

    #start multiprocessing
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Creating consensus Fasta\n')
        #start multiprocessing
        results = pool.starmap_async(cd.call,[[cmd,'/outfile',{reference_sequence_dir:"/ref",outDir:"/outfile"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')

    end = time.time()
    runtime = round(end - start,2)
    print(f'\nSniffles: Finished generating consensus fasta in {runtime} seconds')

    #determine runtime of processes
    end = time.time()
    runtime = round(end - overall_start,2)
    print(f'\nSniffles: Finished generating consensus sequence in {runtime} seconds')
    return out_fasta

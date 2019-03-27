import os
import shlex
import subprocess as sub
import multiprocessing as mp
import time
from sc import procTitle,checkexists
import calldocker as cd

def snpcaller(runCFG,bam_files,threads='1'):
    #set parameters
    outDir =runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])
    outDir = os.path.join(outDir,'snp_calls')
    checkexists(outDir)

    #set reference sequence
    reference_sequence_abspath = os.path.abspath(runCFG['exec']['referenceSequence'])
    reference_sequence_name = os.path.basename(reference_sequence_abspath)
    reference_sequence_dir = os.path.dirname(reference_sequence_abspath)

    #starting time point
    start =  time.time()
    procTitle('SNP Calling')
    if runCFG['snpcalling']['caller'] == 'lofreq':
        caller = 'LoFreq'
    if runCFG['snpcalling']['caller'] == 'varscan':
        caller = 'VarScan'
    print(f'\nSniffles: Started calling SNPs using {caller}')

    cmds = []
    for path in bam_files:
        full_path = os.path.abspath(path)
        file_name = os.path.basename(full_path)
        path = os.path.dirname(full_path)
        id = file_name.split(".")[0]
        #using lofreq
        if caller == 'LoFreq':
            #call snps
            cmd1 = f'bash -c \'lofreq call -f /ref/{reference_sequence_name} -o {id}_lofreq.vcf /infile/{id}.bam && '

            #filter snps
            min_cov = runCFG['snpcalling']['minCoverage']
            snp_qual_threshold = runCFG['snpcalling']['snpQualityThreshold']
            snp_frequency = runCFG['snpcalling']['snpFrequency']
            cmd2 = f'lofreq filter --cov-min {min_cov} --snvqual-thresh {snp_qual_threshold} --af-min {snp_frequency} -i {id}_lofreq.vcf -o {id}_snps.vcf \''

            #add commands to list for multiprocessing
            cmds.append(cmd1+cmd2)

        #using varscan
        if caller == 'VarScan':
            #generate mpileup
            cmd1 = f'bash -c \'samtools mpileup -d 1000000 /infile/{id}.bam -f /ref/{reference_sequence_name} > {id}_snps.pileup &&'

            #call snps
            snp_frequency=runCFG['snpcalling']['snpFrequency']
            min_cov=runCFG['snpcalling']['minCoverage']
            snp_qual_threshold=runCFG['snpcalling']['snpQualityThreshold']

            cmd2 = f'java -jar /tools/varscan.jar mpileup2snp {id}_snps.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1 > {id}_snps.vcf\''

            #add commands to list for multiprocessing
            cmds.append(cmd1+cmd2)

        #command list for annotating aa changes
        if runCFG['exec']['annotateAAChanges']:
            pass
            #TODO add annotater for annotating aa changes
            #reference = reference.split('.')[0]
            #cmd = f'java -jar {libPath}/snpeff/snpEff.jar {reference}.fasta {id}_snps.vcf'
            #outfile = f'{id}_snps_annotated.vcf'
            #annotate_cmds.append([cmd,outfile])

    #execute multiprocessing commands
    #start multiprocessing
    pool = mp.Pool(processes=1)

    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Calling SNPs\n')
        results = pool.starmap_async(cd.call,[[cmd,'/outfile',{reference_sequence_dir:"/ref",path:"/infile",outDir:"/outfile"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')

    #get end time
    end = time.time()
    #get total runtime
    runtime = end - start
    print(f'\nSniffles: Finished calling snps in {runtime} seconds')

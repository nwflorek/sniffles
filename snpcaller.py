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
    print(f'\nSniffles: Started calling SNPs')

    bams = []
    sample_list = []
    for path in bam_files:
        full_path = os.path.abspath(path)
        file_name = os.path.basename(full_path)
        path = os.path.dirname(full_path)
        id = file_name.split(".")[0]
        sample_list.append(id)
        bams.append('/infile/'+file_name)

    #generate mpileup
    cmd1 = 'bash -c \'samtools mpileup -R -d 1000000 {bams} -f /ref/{reference_sequence_name} > all.mpileup &&'.format(bams=' '.join(bams),reference_sequence_name=reference_sequence_name)

    #call snps
    snp_frequency=runCFG['snpcalling']['snpFrequency']
    min_cov=runCFG['snpcalling']['minCoverage']
    snp_qual_threshold=runCFG['snpcalling']['snpQualityThreshold']

    cmd2 = 'java -jar /tools/varscan.jar mpileup2cns all.mpileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1 --variants --vcf-sample-list <(echo -e "{samples}") > all_snps.vcf\''.format(min_cov=min_cov,snp_qual_threshold=snp_qual_threshold,snp_frequency=snp_frequency,samples='\n'.join(sample_list))

    #add commands to list for multiprocessing
    cmd =cmd1+cmd2

    #future code block for annotating aa changes
    if runCFG['exec']['annotateAAChanges']:
        pass
        #TODO add annotater for annotating aa changes

    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Calling SNPs\n')
        results = cd.call(cmd,'/outfile',{reference_sequence_dir:"/ref",path:"/infile",outDir:"/outfile"})
        outlog.write('-----------\n')
        outlog.write(results)
        #denote end of logs
        outlog.write('***********\n')

    #get end time
    end = time.time()
    #get total runtime
    runtime = round(end - start,2)
    print(f'\nSniffles: Finished calling snps in {runtime} seconds')

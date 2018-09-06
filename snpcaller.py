import os
import shlex
import subprocess as sub
import multiprocessing as mp
import time
from sniffProc import proc,init
from sc import procTitle

def snpcaller(readData,runCFG,threads='1',ids=''):
    #easy const access
    libPath = runCFG['libPath']
    outDir =runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])

    #use id if provided otherwise get list
    if not ids:
        ids = readData.idList

    #initalize command lists
    lofreq_cmds = []
    varscan_cmds = []
    annotate_cmds = []

    #starting time point
    start =  time.time()
    procTitle('SNP Calling')
    if runCFG['snpcalling']['lofreq']:
        caller = 'LoFreq'
    else:
        caller = 'VarScan'
    print(f'\nSniffles: Started calling SNPs using {caller}')

    for id in ids:
        #determine samfile that will be used
        if runCFG['exec']['removeDupReads']:
            samfile = f'{id}_nodups.sam'
        elif runCFG['exec']['normalizeCoverage']:
            samfile = f'{id}_remapped.sam'
        else:
            samfile = f'{id}.sam'

        #convert sam to sorted bam file
        cmd01 = f'{libPath}/bin/samtools view -b {samfile}'
        cmd01 = shlex.split(cmd01)
        cmd02 = f'{libPath}/bin/samtools sort'
        cmd02 = shlex.split(cmd02)
        with open(runCFG['exec']['outdir']+'/'+f'{id}_snps.bam','w') as outbam:
            c1 = sub.Popen(cmd01,stdout=sub.PIPE,cwd=runCFG['exec']['outdir'])
            c2 = sub.Popen(cmd02,stdin=c1.stdout,stdout=outbam,cwd=runCFG['exec']['outdir'])
            c2.wait()

        #call snps and annotate
        reference = runCFG['exec']['referenceSequence']

        #using lofreq
        if runCFG['snpcalling']['lofreq']:
            #call snps
            cmd01 = f'{libPath}/bin/lofreq call -f {reference} -o {id}_lofreq.vcf {id}.bam'
            cmd01 = shlex.split(cmd)

            #filter snps
            min_cov = runCFG['snpcalling']['minCoverage']
            snp_qual_threshold = runCFG['snpcalling']['snpQualityThreshold']
            snp_frequency = runCFG['snpcalling']['snpFrequency']
            cmd02 = f'{libPath}/bin/lofreq filter --cov-min {min_cov} --snvqual-thresh {snp_qual_threshold} --af-min {snp_frequency} -i {id}_lofreq.vcf -o {id}_snps.vcf'
            cmd02 = shlex.split(cmd)

            #add commands to list for multiprocessing
            lofreq_cmds.append([cmd01,cmd02])

        #using varscan
        if runCFG['snpcalling']['varscan']:
            #generate mpileup
            cmd01 = f'{libPath}/bin/samtools mpileup -d 1000000 {id}.bam -f {reference}'
            outfile01 = f'{id}_snps.pileup'

            #call snps
            snp_frequency=runCFG['snpcalling']['snpFrequency']
            min_cov=runCFG['snpcalling']['minCoverage']
            snp_qual_threshold=runCFG['snpcalling']['snpQualityThreshold']
            cmd02 = f'java -jar {libPath}/varscan/VarScan.v2.3.9.jar mpileup2snp {id}_snps.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1'
            outfile02 = f'{id}_snps.vcf'

            #add commands to list for multiprocessing
            varscan_cmds.append([cmd01,outfile01,cmd02,outfile02])

        #command list for annotating aa changes
        if runCFG['exec']['annotateAAChanges']:
            #TODO add annotater for annotating aa changes
            #reference = reference.split('.')[0]
            #cmd = f'java -jar {libPath}/snpeff/snpEff.jar {reference}.fasta {id}_snps.vcf'
            #outfile = f'{id}_snps_annotated.vcf'
            #annotate_cmds.append([cmd,outfile])

    #execute multiprocessing commands
    #start multiprocessing
    lock = mp.Lock()
    pool = mp.Pool(processes=1,initializer=init,initargs=(lock,))

    if lofreq_cmds:
        pool.starmap(proc, [[runCFG,i[0]] for i in lofreq_cmds])
        pool.starmap(proc, [[runCFG,i[1]] for i in lofreq_cmds])
    if varscan_cmds:
        pool.starmap(proc, [[runCFG,i[0],'',i[1]] for i in varscan_cmds])
        pool.starmap(proc, [[runCFG,i[2],'',i[3]] for i in varscan_cmds])
    if annotate_cmds:
        pool.starmap(proc, [[runCFG,i[0],'',i[1]] for i in annotate_cmds])

    #get end time
    end = time.time()
    #get total runtime
    runtime = end - start
    print(f'\nSniffles finished calling snps in {runtime} seconds')

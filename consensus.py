import os
import shlex
import subprocess as sub
import mapping

def consensus(readData,runCFG,threads,ids=''):
    #easy const access
    libPath = runCFG['libPath']
    outDir =runCFG['exec']['outdir']

    #if no id supplied get list of all ids
    if not ids:
        ids = readData.idList

    #add consensus extraction tracking
    if 'consensus' not in readData.data:
          readData.data['consensus'] = []

    #if no consensus dir exists create one
    try:
        os.mkdir(os.path.join(outDir,'consensus'))
    except:
        pass

    #set reference sequence and begin
    reference = runCFG['exec']['referenceSequence']
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
        with open(outDir+'/'+f'{id}.bam','w') as outbam:
            c1 = sub.Popen(cmd01,stdout=sub.PIPE,cwd=outDir)
            c2 = sub.Popen(cmd02,stdin=c1.stdout,stdout=outbam,cwd=outDir)
            c2.wait()

        #make multiway pileup using samtools
        print(f'generating consensus on {id} using varscan')
        cmd = f'{libPath}/bin/samtools mpileup -d 1000000 {id}.bam -f {reference}'
        cmd = shlex.split(cmd)
        with open(outDir+'/'+f'{id}.pileup','w') as mpileup:
            sub.Popen(cmd,cwd=outDir,stdout=mpileup).wait()

        #run varscan mpileup2cns to generate vcf with consensus information
        minCov = runCFG['snpcalling']['minCoverage']
        quality = runCFG['snpcalling']['snpQualityThreshold']
        freq = runCFG['snpcalling']['consensusFrequency']
        cmd = f'java -jar {libPath}/varscan/VarScan.v2.3.9.jar mpileup2cns {id}.pileup --min-coverage {minCov} --min-avg-qual {quality} --min-var-freq {freq} --strand-filter 1 --output-vcf 1'
        cmd = shlex.split(cmd)
        with open(outDir+'/consensus/'+f'{id}.vcf','w') as vcfout:
            sub.Popen(cmd,cwd=outDir,stdout=vcfout).wait()

        #compress vcf file with bgzip and index with tabix
        cmd = f'{libPath}/bin/bgzip consensus/{id}.vcf'
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=outDir).wait()
        cmd = f'{libPath}/bin/tabix consensus/{id}.vcf.gz'
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=outDir).wait()

        #use bcftools to get consensus fasta
        cmd = f'{libPath}/bin/bcftools consensus -f {reference} consensus/{id}.vcf.gz'
        cmd = shlex.split(cmd)
        with open(outDir+'/consensus/'+f'{id}.fasta','w') as fastaout:
            sub.Popen(cmd,cwd=outDir,stdout=fastaout).wait()

        #add id as finished on tracking
        readData.data['consensus'].append(id)

    #if mapping reads to consensus is specified run mapping
    if runCFG['exec']['mapToConsensus']:
        mapping(readData,runCFG,threads,ids)

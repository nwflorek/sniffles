import os;
import shlex
import subprocess as sub

def consensus(libpath,runCFG,threads,ids):
    os.mkdir(runCFG['exec']['outdir']+'/'+'consensus')
    reference = runCFG['exec']['referenceSequence']
    for id in ids:
        #determine samfile that will be used
        if runCFG['exec']['removeDupReads']:
            samfile = '{id}_nodups.sam'.format(id=id)
        elif runCFG['exec']['normalizeCoverage']:
            samfile = '{id}_remapped.sam'.format(id=id)
        else:
            samfile = '{id}.sam'.format(id=id)

        #convert sam to sorted bam file
        cmd01 = '{libpath}/bin/samtools view -b {samfile}'.format(libpath=libpath,id=id,samfile=samfile)
        cmd01 = shlex.split(cmd01)
        cmd02 = '{libpath}/bin/samtools sort'.format(libpath=libpath)
        cmd02 = shlex.split(cmd02)
        with open(runCFG['exec']['outdir']+'/'+'{id}.bam'.format(id=id),'w') as outbam:
            c1 = sub.Popen(cmd01,stdout=sub.PIPE,cwd=runCFG['exec']['outdir'])
            c2 = sub.Popen(cmd02,stdin=c1.stdout,stdout=outbam,cwd=runCFG['exec']['outdir'])
            c2.wait()

        print('generating consensus on {id} using varscan'.format(id=id))
        cmd = "{libpath}/bin/samtools mpileup -d 1000000 {id}.bam -f {reference_sequence}".format(libpath=libpath, id=id, reference_sequence=reference)
        cmd = shlex.split(cmd)
        with open(runCFG['exec']['outdir']+'/'+'{id}.pileup'.format(id=id),'w') as mpileup:
            sub.Popen(cmd,cwd=runCFG['exec']['outdir'],stdout=mpileup).wait()

        cmd = "java -jar {libpath}/varscan/VarScan.v2.3.9.jar mpileup2cns {id}.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq 0.5 --strand-filter 1 --output-vcf 1".format(libpath=libpath, id=id, snp_frequency=runCFG['snpcalling']['snpFrequency'],min_cov=runCFG['snpcalling']['minCoverage'], snp_qual_threshold=runCFG['snpcalling']['snpQualityThreshold'])
        cmd = shlex.split(cmd)
        with open(runCFG['exec']['outdir']+'/consensus/'+'{id}.vcf'.format(id=id),'w') as vcfout:
            sub.Popen(cmd,cwd=runCFG['exec']['outdir'],stdout=vcfout).wait()

        #compress vcf file with bgzip and index with tabix
        cmd = '{libpath}/bin/bgzip consensus/{id}.vcf'.format(libpath=libpath,id=id)
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()
        cmd = '{libpath}/bin/tabix consensus/{id}.vcf.gz'.format(libpath=libpath,id=id)
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()

        #use bcftools to get consensus fasta
        cmd = "{libpath}/bin/bcftools consensus -f {reference} consensus/{id}.vcf.gz".format(libpath=libpath,reference=reference,id=id)
        cmd = shlex.split(cmd)
        with open(runCFG['exec']['outdir']+'/consensus/'+'{id}.fasta'.format(id=id),'w') as fastaout:
            sub.Popen(cmd,cwd=runCFG['exec']['outdir'],stdout=fastaout).wait()

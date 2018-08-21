import subprocess as sub
import os
import shlex

def snpcaller(libpath,runCFG,threads,ids):
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

        #call snps and annotate
        reference = runCFG['exec']['referenceSequence']
        #using lofreq
        if runCFG['snpcalling']['lofreq']:
            print('now calling SNPs on {id} using lofreq'.format(id=id))
            cmd = "{libpath}/bin/lofreq call -f {reference_sequence} -o {id}_lofreq.vcf {id}.bam".format(libpath=libpath,id=id,reference_sequence=reference)
            cmd = shlex.split(cmd)
            sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()
            cmd = "{libpath}/bin/lofreq filter --cov-min {min_cov} --snvqual-thresh {snp_qual_threshold} --af-min {snp_frequency} -i {id}_lofreq.vcf -o {id}_lofreq_filtered.vcf".format(libpath=libpath,id=id, min_cov=runCFG['snpcalling']['minCoverage'], snp_qual_threshold=runCFG['snpcalling']['snpQualityThreshold'], snp_frequency=runCFG['snpcalling']['snpFrequency'])
            cmd = shlex.split(cmd)
            sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()

            if runCFG['exec']['annotateAAChanges']:
                cmd = "java -jar {libpath}/snpeff/snpEff.jar {reference_sequence} {id}_lofreq_filtered.vcf".format(libpath=libpath,id=id,reference_sequence=reference)
                cmd = shlex.split(cmd)
                with open(runCFG['exec']['outdir']+'/'+'{id}_lofreq_annotated.vcf'.format(id=id),'w') as outvcf:
                    sub.Popen(cmd,cwd=runCFG['exec']['outdir'],stdout=outvcf).wait()

        #using varscan
        if runCFG['snpcalling']['varscan']:
            print('now calling SNPs on {id} using varscan'.format(id=id))
            cmd = "{libpath}/bin/samtools mpileup -d 1000000 {id}.bam -f {reference_sequence}".format(libpath=libpath, id=id, reference_sequence=reference)
            cmd = shlex.split(cmd)
            with open(runCFG['exec']['outdir']+'/'+'{id}.pileup'.format(id=id),'w') as mpileup:
                sub.Popen(cmd,cwd=runCFG['exec']['outdir'],stdout=mpileup).wait()

            cmd = "java -jar {libpath}/varscan/VarScan.v2.3.9.jar mpileup2snp {id}.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1".format(libpath=libpath, id=id, snp_frequency=runCFG['snpcalling']['snpFrequency'],min_cov=runCFG['snpcalling']['minCoverage'], snp_qual_threshold=runCFG['snpcalling']['snpQualityThreshold'])
            cmd = shlex.split(cmd)
            with open(runCFG['exec']['outdir']+'/'+'{id}_varscan.vcf'.format(id=id),'w') as vcfout:
                sub.Popen(cmd,cwd=runCFG['exec']['outdir'],stdout=vcfout).wait()

            if runCFG['exec']['annotateAAChanges']:
                cmd = "java -jar {libpath}/snpeff/snpEff.jar {reference_sequence} {id}_varscan.vcf".format(libpath=libpath,id=id,reference_sequence=reference.split('.')[0])
                cmd = shlex.split(cmd)
                with open(runCFG['exec']['outdir']+'/'+'{id}_varscan_annotated.vcf'.format(id=id),'w') as outvcf:
                    sub.Popen(cmd,cwd=runCFG['exec']['outdir'],stdout=outvcf).wait()

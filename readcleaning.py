import subprocess as sub
import os
import shlex
from mapping import mapping

def removeDuplicates(libpath,runCFG,threads,ids):
    for id in ids:
        #sort samfile
        cmd01 = '{libpath}/bin/samtools view -b {id}.sam'.format(libpath=libpath,id=id)
        cmd01 = shlex.split(cmd01)
        cmd02 = '{libpath}/bin/samtools sort'.format(libpath=libpath)
        cmd02 = shlex.split(cmd02)
        cmd03 = '{libpath}/bin/samtools view -h'.format(libpath=libpath)
        cmd03 = shlex.split(cmd03)
        with open(runCFG['exec']['outdir']+'/'+'{id}_sorted.sam'.format(id=id),'w') as outsam:
            c1 = sub.Popen(cmd01,stdout=sub.PIPE,cwd=runCFG['exec']['outdir'])
            c2 = sub.Popen(cmd02,stdin=c1.stdout,stdout=sub.PIPE,cwd=runCFG['exec']['outdir'])
            c3 = sub.Popen(cmd03,stdin=c2.stdout,stdout=outsam,cwd=runCFG['exec']['outdir'])
            c3.wait()

        #remove duplicate reads
        cmd = 'java -jar {libpath}/picard/picard.jar MarkDuplicates I={id}_sorted.sam O={id}_nodups.sam REMOVE_DUPLICATES=true M={id}.removeDupMetrics.txt'.format(libpath=libpath,id=id)
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()
        os.remove(runCFG['exec']['outdir']+'/'+'{id}_sorted.sam'.format(id=id))
        os.remove(runCFG['exec']['outdir']+'/'+'{id}.sam'.format(id=id))

def normCoverage(libpath,runCFG,threads,ids):

    for id in ids:
        #determine which samfile to use if duplicates have been removed
        if runCFG['exec']['removeDupReads']:
            samfile = '{id}_nodups.sam'.format(id=id)
        else:
            samfile = '{id}.sam'.format(id=id)

        #get reads from samfile
        cmd = '{libpath}/bbmap/reformat.sh in={samfile} out={id}_adjusted.fastq'.format(libpath=libpath,samfile=samfile,id=id)
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()

        #run bbnorm
        cmd = '{libpath}/bbmap/bbnorm.sh in={id}_adjusted.fastq out={id}_normalized.fastq target={depth}'.format(libpath=libpath,id=id,depth=runCFG['exec']['coverageNormDepth'])
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()
        os.remove('{path}/{id}_adjusted.fastq'.format(path=runCFG['exec']['outdir'],id=id))

    mapping(libpath,runCFG,threads,ids,True)

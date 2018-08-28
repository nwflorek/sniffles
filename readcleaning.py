import subprocess as sub
import os
import shlex
from mapping import mapping

def removeDuplicates(readData,runCFG,threads='1',ids=''):
    #easy const access
    libPath = runCFG['libPath']
    outDir =runCFG['exec']['outdir']

    #if no id supplied get list of all ids
    if not ids:
        ids = readData.idList

    for id in ids:
        #sort samfile
        cmd01 = f'{libPath}/bin/samtools view -b {id}.sam'
        cmd01 = shlex.split(cmd01)
        cmd02 = f'{libpath}/bin/samtools sort'
        cmd02 = shlex.split(cmd02)
        cmd03 = f'{libpath}/bin/samtools view -h'
        cmd03 = shlex.split(cmd03)
        with open(outDir+'/'+'{id}_sorted.sam'.format(id=id),'w') as outsam:
            c1 = sub.Popen(cmd01,stdout=sub.PIPE,cwd=outDir)
            c2 = sub.Popen(cmd02,stdin=c1.stdout,stdout=sub.PIPE,cwd=outDir)
            c3 = sub.Popen(cmd03,stdin=c2.stdout,stdout=outsam,cwd=outDir)
            c3.wait()

        #remove duplicate reads
        cmd = f'java -jar {libPath}/picard/picard.jar MarkDuplicates I={id}_sorted.sam O={id}_nodups.sam REMOVE_DUPLICATES=true M={id}.removeDupMetrics.txt'
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=outDir).wait()
        os.remove(outDir+'/'+'{id}_sorted.sam'.format(id=id))
        os.remove(outDir+'/'+'{id}.sam'.format(id=id))

        #add id to finished list
        readData.addData('rmDuplicates',id,f'{id}_nodups.sam')

def normCoverage(readData,runCFG,threads='1',ids=''):
    #easy const access
    libPath = runCFG['libPath']
    outDir =runCFG['exec']['outdir']

    #if no id supplied get list of all ids
    if not ids:
        ids = readData.idList

    for id in ids:
        #determine which samfile to use if duplicates have been removed
        if id in readData.data['rmDuplicates']:
            samfile = f'{id}_nodups.sam'
        else:
            samfile = f'{id}.sam'

        #get reads from samfile
        cmd = f'{libPath}/bbmap/reformat.sh in={samfile} out={id}_adjusted.fastq'
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=outDir).wait()

        #run bbnorm
        cov = runCFG['exec']['coverageNormDepth']
        cmd = f'{libPath}/bbmap/bbnorm.sh in={id}_adjusted.fastq out={id}_normalized.fastq target={cov}'
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=outDir).wait()
        os.remove(f'{outDir}/{id}_adjusted.fastq')

        #add to tracker
        readData.addData('normalized',id,f'{outDir}/{id}_normalized.fastq')

    mapping(readData,runCFG,threads,ids)

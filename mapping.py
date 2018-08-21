import os;
import shlex
import subprocess as sub

def mapping(libPath,runCFG,threads,ids,normalized=False):
    refPath = os.path.abspath(runCFG['exec']['referenceSequence'])
    refName = os.path.basename(refPath)

    #index reference
    cmd = '{path}/bin/bowtie2-build {reference_sequence} {reference_sequence_name}'.format(path=libPath,reference_sequence=refPath,reference_sequence_name=refName)
    cmd = shlex.split(cmd)
    sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()

    #map reads
    fwdReads = []
    revReads = []
    if normalized:
        for id in ids:
            reads = '{id}_normalized.fastq'.format(id=id)
            cmd = "{path}/bin/bowtie2 -x {reference_sequence} --interleaved {reads} -S {id}_remapped.sam -p {threads} --local".format(path=libPath,id=id,reads=reads,threads=threads,reference_sequence=refName)
            cmd = shlex.split(cmd)
            sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()
    else:
        for id in ids:
            fwdRead = 'trimmed/{id}_trimmed_1P.fastq.gz'.format(id=id)
            revRead = 'trimmed/{id}_trimmed_2P.fastq.gz'.format(id=id)
            cmd = "{path}/bin/bowtie2 -x {reference_sequence} -1 {fwdRead} -2 {revRead} -S {id}.sam -p {threads} --local".format(path=libPath,id=id,fwdRead=fwdRead,revRead=revRead,threads=threads,reference_sequence=refName)
            cmd = shlex.split(cmd)
            sub.Popen(cmd,cwd=runCFG['exec']['outdir']).wait()

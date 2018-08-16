import os
import shlex
import subprocess as sub
def trimmomatic(libPath,runCFG,id,read1,read2):
    #parameters
    minlength = runCFG['trimmomatic']['minlength']
    windowsize = runCFG['trimmomatic']['windowSize']
    qscore = runCFG['trimmomatic']['qscore']
    adapterpath = runCFG['trimmomatic']['adaptersPath']

    #main command
    cmd = 'java -jar {libPath}/trimmomatic/trimmomatic-0.36.jar '.format(libPath=libPath)

    #determine args
    if runCFG['trimmomatic']['removeAdapters']:
        if runCFG['trimmomatic']['paired']:
            args = 'PE {read1} {read2} -baseout {id}_trimmed.fastq.gz ILLUMINACLIP:{adapters_fasta}:1:30:10 SLIDINGWINDOW:{window}:{qscore} MINLEN:{minlen}'.format(id=id,read1=read1,read2=read2,window=windowsize,qscore=qscore,minlen=minlength,adapters_fasta=adapterpath)
        else:
            args = 'SE {read} -baseout {id}_trimmed.fastq.gz ILLUMINACLIP:{adapters_fasta}:1:30:10 SLIDINGWINDOW:{window}:{qscore} MINLEN:{minlen}'.format(id=id,read=read1,window=windowsize,qscore=qscore,minlen=minlength,adapters_fasta=adapterpath)
    else:
        if runCFG['trimmomatic']['paired']:
            args = 'PE {read1} {read2} -baseout {id}_trimmed.fastq.gz SLIDINGWINDOW:{window}:{qscore} MINLEN:{minlen}'.format(id=id,read1=read1,read2=read2,window=windowsize,qscore=qscore,minlen=minlength)
        else:
            args = 'SE {read} -baseout {id}_trimmed.fastq.gz SLIDINGWINDOW:{window}:{qscore} MINLEN:{minlen}'.format(id=id,read=read1,window=windowsize,qscore=qscore,minlen=minlength)
    cmd = shlex.split(cmd+args)
    try:
        os.mkdir(os.path.join(runCFG['exec']['outdir'],'trimmed'))
    except:
        pass
    sub.Popen(cmd,cwd=os.path.join(runCFG['exec']['outdir'],'trimmed')).wait()

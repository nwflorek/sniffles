import os
import shlex
import subprocess as sub

#define process for multiprocessing
def proc(runCFG,cmd,d='',outfile='',shell=False):
    #set parameters
    outDir = runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])
    
    if shell:
        p = sub.Popen(cmd,cwd=os.path.join(outDir,d),shell=True,stdout=sub.PIPE,stderr=sub.PIPE)
        out,err = p.communicate()
    else:
        cmd = shlex.split(cmd)

        #start process
        if outfile:
            with open(os.path.join(outDir,d,outfile),'w') as out:
                p = sub.Popen(cmd,cwd=os.path.join(outDir,d),stdout=out,stderr=sub.PIPE)
                err = p.stderr.read()
            out = ''

        else:
            p = sub.Popen(cmd,cwd=os.path.join(outDir,d),stdout=sub.PIPE,stderr=sub.PIPE)
            out,err = p.communicate()

    #open logfile and store log information
    lock.acquire()
    with open(logfile,'a') as outlog:
        if out:
            out = out.decode('utf-8')
            for line in out:
                outlog.write(line)
        if err:
            err = err.decode('utf-8')
            for line in err:
                outlog.write(line)
    lock.release()

#make a lock global for child workers
def init(l):
    global lock
    lock = l

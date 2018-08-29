import os
import shlex
import subprocess as sub

#define process for multiprocessing
def proc(runCFG,cmd,d=''):
    #set parameters
    outDir = runCFG['exec']['outdir']
    logfile = os.path.join(outDir,runCFG['exec']['logfile'])
    cmd = shlex.split(cmd)

    #start process
    p = sub.Popen(cmd,cwd=os.path.join(outDir,d),stdout=sub.PIPE,stderr=sub.PIPE)
    out,err = p.communicate()

    #open logfile and store log information
    lock.acquire()
    with open(logfile,'a') as outlog:
        out = out.decode('utf-8')
        err = err.decode('utf-8')
        if out:
            for line in out:
                outlog.write(line)
        if err:
            for line in err:
                outlog.write(line)
    lock.release()

#make a lock global for child workers
def init(l):
    global lock
    lock = l

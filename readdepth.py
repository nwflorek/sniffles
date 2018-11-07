import subprocess as sub
import sys
import os
import shlex

def average_depth(path):
    if '.bam' in path:
        cmd = f'samtools depth -a {path} | awk \'{{sum+=$3}} END {{ print "Average:",sum/NR}}\''
    else:
        return 'error: not a bam file'

    p = sub.Popen(cmd,shell=True,stdout=sub.PIPE,stderr=sub.PIPE)
    out,err = p.communicate()
    try:
        depth = out.decode('utf-8').strip().split(':')[1][1:]
    except IndexError:
        o = out.decode('utf-8')
        return f'error: {o}'
    return depth

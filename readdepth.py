import subprocess as sub
import sys
import os
import shlex
import calldocker as cd

def average_depth(file):
    abs_path = os.path.abspath(file)
    path = os.path.dirname(abs_path)
    bam = os.path.basename(abs_path)

    if '.bam' in bam:
        cmd = f'bash -c "samtools depth -a {bam} | awk \'{{sum+=$3}} END {{ print sum/NR}}\'"'
    else:
        return f'error: {bam} is not a bam file'
    out = cd.call(cmd,'/data',{path:"/data"})
    return out

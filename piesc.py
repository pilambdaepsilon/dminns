import subprocess, shlex
import sys
import time

def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print proc_stdout
m = 1
mstr = str(m)
subprocess_cmd('if [ ! -d \"BOSONS/DmNoFiles/%sGeV\" ]; then echo "CREATING DIRECTORY BOSONS/DmNoFiles/%sGeV..."; mkdir \"BOSONS/DmNoFiles/%sGeV\"; fi' %(mstr,mstr,mstr))

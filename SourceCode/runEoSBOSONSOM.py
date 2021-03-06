import subprocess
import numpy as np
import sys
import os
import time

to = time.time()
def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print proc_stdout

startit = float(sys.argv[1])
endit = float(sys.argv[2])

startsigchinpow = sys.argv[3]
endsigchinpow = sys.argv[4]

startmass = sys.argv[6]
endmass = sys.argv[7]

sigchi2pow = int(sys.argv[5])

STARNAME = sys.argv[8]
MODELNAME = sys.argv[9]

massinc = sys.argv[10]

powstrchi2 = str(sigchi2pow)
#m = float(startmass)
#while m <= float(endmass):
m = float(startmass)
if(m>=1):
    m = int(startmass)
while m <= float(endmass):
    mstr = str(m)
    subprocess_cmd('if [ ! -d \"BOSONS/DmNoFiles/%sGeV\" ]; then echo "CREATING DIRECTORY BOSONS/DmNoFiles/%sGeV..."; mkdir \"BOSONS/DmNoFiles/%sGeV\"; fi' %(mstr,mstr,mstr))
    subprocess_cmd('if [ ! -d \"BOSONS/EoSFiles/%sGeV\" ]; then echo "CREATING DIRECTORY BOSONS/EoSFiles/%sGeV..."; mkdir \"BOSONS/EoSFiles/%sGeV\"; fi' %(mstr,mstr,mstr))
    subprocess_cmd('if [ ! -d \"BOSONS/MRFiles/%sGeV\" ]; then echo "CREATING DIRECTORY BOSONS/MRFiles/%sGeV..."; mkdir \"BOSONS/MRFiles/%sGeV\"; fi' %(mstr,mstr,mstr))
#===========================================================
#================== RUN No. Dens. FILES ====================
#===========================================================
    for SIGCHIN in range (int(startsigchinpow), int(endsigchinpow)):
        SIGCHINstr = str(SIGCHIN)
        print("sigchin: %s"%SIGCHINstr)
        print("sigchi2: %s"%powstrchi2)
        print("mass: %s"%mstr)
        #RUN NUMBER DENSITY FILES
#===========================================================
#================== RUN EoS FILES ==========================
#===========================================================
        d = startit
        while(d < endit):
            sys.stdout.write("\r ./EoSDMBOM %.1f %s %s %s %s %s" %(d, SIGCHINstr, powstrchi2, mstr, STARNAME, MODELNAME))
            sys.stdout.flush()
            dstr = str(d)
            proc = subprocess.check_output(["./EoSDMBOM", "%.1f"%d, "%s"%SIGCHINstr, "%s"%powstrchi2, "%s"%mstr, "%s"%STARNAME, "%s"%MODELNAME], stderr=subprocess.STDOUT) 
            d += 0.1
        print("\nEoS done...\n")

    if m>1:
        m += int(massinc)
    else:
        m += float(massinc)
#        m += 0.1


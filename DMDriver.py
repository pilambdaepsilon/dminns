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

sigchi2pow = int(sys.argv[5])

startmass = sys.argv[6]
endmass = sys.argv[7]

STARNAME = sys.argv[8]
MODELNAME = sys.argv[9]

massinc = sys.argv[10]

PARTICLE = sys.argv[11]

powstrchi2 = str(sigchi2pow)
#m = float(startmass)
#while m <= float(endmass):
m = float(startmass)
if(m>=1):
    m = int(startmass)
while m <= float(endmass):
    mstr = str(m)
    subprocess_cmd('if [ ! -d \"%s/DmNoFiles/%sGeV\" ]; then echo "CREATING DIRECTORY %s/DmNoFiles/%sGeV..."; mkdir \"%s/DmNoFiles/%sGeV\"; fi' %(PARTICLE,mstr,PARTICLE,mstr,PARTICLE,mstr))
    subprocess_cmd('if [ ! -d \"%s/EoSFiles/%sGeV\" ]; then echo "CREATING DIRECTORY %s/EoSFiles/%sGeV..."; mkdir \"%s/EoSFiles/%sGeV\"; fi' %(PARTICLE,mstr,PARTICLE,mstr,PARTICLE,mstr))
    subprocess_cmd('if [ ! -d \"%s/MRFiles/%sGeV\" ]; then echo "CREATING DIRECTORY %s/MRFiles/%sGeV..."; mkdir \"%s/MRFiles/%sGeV\"; fi' %(PARTICLE,mstr,PARTICLE,mstr,PARTICLE,mstr))
#===========================================================
#================== RUN No. Dens. FILES ====================
#===========================================================
    for SIGCHIN in range (int(startsigchinpow), int(endsigchinpow)+1):
        SIGCHINstr = str(SIGCHIN)
        print("sigchin: %s"%SIGCHINstr)
        print("sigchi2: %s"%powstrchi2)
        print("mass: %s"%mstr)
        #RUN NUMBER DENSITY FILES
        d = startit
        while(d <= endit):
            sys.stdout.write("\r ./exe/DMNumberDensity %.1f %s %s %s %s %s %s" %(d, SIGCHINstr, powstrchi2, mstr, STARNAME, MODELNAME, PARTICLE))
            sys.stdout.flush()
            dstr = str(d)
            proc = subprocess.check_output(["./exe/DMNumberDensity", "%.1f"%d, "%s"%SIGCHINstr, "%s"%powstrchi2, "%s"%mstr, "%s"%STARNAME, "%s"%MODELNAME, "%s"%PARTICLE], stderr=subprocess.STDOUT) 
            d += 0.1
        print("\nNumber Density done...\n")
#===========================================================
#================== RUN EoS FILES ==========================
#===========================================================
        d = startit
        while(d <= endit):
            sys.stdout.write("\r ./exe/DMEquationOfState %.1f %s %s %s %s %s %s" %(d, SIGCHINstr, powstrchi2, mstr, STARNAME, MODELNAME, PARTICLE))
            sys.stdout.flush()
            dstr = str(d)
            proc = subprocess.check_output(["./exe/DMEquationOfState", "%.1f"%d, "%s"%SIGCHINstr, "%s"%powstrchi2, "%s"%mstr, "%s"%STARNAME, "%s"%MODELNAME, "%s"%PARTICLE], stderr=subprocess.STDOUT) 
            d += 0.1
        print("\nEoS done...\n")

#===========================================================
#================== RUN M-R FILES ==========================
#===========================================================
        d = startit
        while(d <= endit):
            dt = (time.time() - to)/60
            sys.stdout.write("\r ./exe/DMTOV %.1f %s %s %s %s %s %s... TIME: %.2f min." %(d, SIGCHINstr, powstrchi2, mstr, STARNAME, MODELNAME, PARTICLE, dt))
            sys.stdout.flush()
            dstr = str(d)
            #proc = subprocess.check_output(["./TOVBOM", "%.1f"%d, "%s"%SIGCHINstr, "%s"%powstrchi2, "%s"%mstr, "%s"%STARNAME, "%s"%MODELNAME], stderr=subprocess.STDOUT) 
            retval = os.system("./exe/DMTOV %.1f %s %s %s %s %s %s" %(d, SIGCHINstr, powstrchi2, mstr, STARNAME, MODELNAME, PARTICLE))
            d += 0.1
        print("\nMR done...\n")
    if m>1:
        m += int(massinc)
    else:
        m += float(massinc)
#        m += 0.1

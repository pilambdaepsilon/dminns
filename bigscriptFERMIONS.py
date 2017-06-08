import subprocess
import sys
import time

to = time.time()

startit = int(sys.argv[1])
endit = int(sys.argv[2])

startsigchinpow = sys.argv[3]
endsigchinpow = sys.argv[4]

sigchi2pow = sys.argv[5]

startmass = sys.argv[6]
endmass = sys.argv[7]

powstrchi2 = str(sigchi2pow)

for m in range (int(startmass), int(endmass)):
    mstr = str(m)
    for SIGCHIN in range (int(startsigchinpow), int(endsigchinpow)):
        SIGCHINstr = str(SIGCHIN)
        print("sigchin: %s"%SIGCHINstr)
        print("sigchi2: %s"%powstrchi2)
        print("mass: %s"%mstr)
        #RUN NUMBER DENSITY FILES
        for d in range (startit,endit):
            dstr = str(d)
            proc = subprocess.check_output(["./DmNoF", "%s"%dstr, "%s"%SIGCHINstr, "%s"%powstrchi2, "%s"%mstr], stderr=subprocess.STDOUT) 
        print("\nNumber Density done...\n")
        #RUN EOS FILES
        for d in range (startit,endit):
            dstr = str(d)
            proc = subprocess.check_output(["./EoSDMF", "%s"%dstr, "%s"%SIGCHINstr, "%s"%powstrchi2, "%s"%mstr], stderr=subprocess.STDOUT) 
        print("EoS done...\n")
        #RUN MR FILES
        for d in range (startit,endit):
            dstr = str(d)
            proc = subprocess.check_output(["./TOVF", "%s"%dstr, "%s"%SIGCHINstr, "%s"%powstrchi2, "%s"%mstr], stderr=subprocess.STDOUT) 
            print(d)
            dt = (time.time() - to)/60
            print("time: %f min" %(dt))
        print("MR done...\n")


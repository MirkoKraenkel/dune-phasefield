#! /usr/bin/env python
import subprocess,sys
number=int(sys.argv[1])
nb=str(number)
#subprocess.call(['./phasefield scheme:test&> pmatrix'+nb+'.out'],shell=True)
subprocess.call(['./phasefieldcoupling scheme:test&> pmatrix'+nb+'.out'],shell=True)
subprocess.call(['./phasefieldfd scheme:test&> pfd'+nb+'.out'],shell=True)
subprocess.call(['kdiff3 pmatrix'+nb+'.out pfd'+nb+'.out'],shell=True)

#! /usr/bin/env python
import subprocess,sys
programms=[ 'phasefield', 'phasefieldfd', 'phasefieldfdrho', 'phasefieldfdlambda' , 'phasefieldmfree', 'phasefieldoem']
begin=int(sys.argv[1])
end=int(sys.argv[2])
for p in programms[begin : end ]:
  outfile=p+'_make.out' 
  c = subprocess.call(['make'+' '+p+' &>'+outfile], shell=True)
  print( c )
  if c == 0:
    subprocess.call(['rm '+outfile],shell=True)

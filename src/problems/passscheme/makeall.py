#! /usr/bin/env python
import subprocess , sys
polorder = sys.argv[1]
print(polorder)
progdir='prog_P'+polorder
programms=[ 'phasefield','phasefield_nc', 'phasefield_nctr']
for p in programms[ : ]:
  outfile=p+'_make.out' 
  c = subprocess.call(['make POLORDER='+polorder+' ' +p+' &>'+outfile], shell=True)
  subprocess.call(['mv '+p+' '+p+'_P'+polorder], shell=True)
  print( c )
  if c == 0:
    subprocess.call(['rm '+outfile],shell=True)

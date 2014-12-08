#! /usr/bin/env python
import subprocess,sys
programms=[ 'phasefield', 'phasefieldcoupling','phasefieldfd', 'phasefieldfdrho', 'phasefieldfdlambda' , 'phasefieldmfree', 'phasefieldoem']
prognumber=int(sys.argv[1])
delta=float(sys.argv[2])
runs=int(sys.argv[3])
startLevel=3
timeStep=1e-4
printCount=1000
prog=programms[prognumber]
while runs > 0:
  subprocess.call(['./'+prog+' parameter_script delta:'+str(delta)+' scheme:'+prog+' '+str(delta)+' timeStep:'+str(timeStep)
  +' startLevel:'+str(startLevel)+' printCount:'+str(printCount)],shell=True)
  delta=0.5*delta
  startLevel=startLevel+1
  timeStep=timeStep*0.5
  printCount=printCount*2
  runs=runs-1;
  

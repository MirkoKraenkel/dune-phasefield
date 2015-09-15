#! /usr/bin/env python
import os,subprocess, sys
namelist = [' helmholtz',' pressure', ' a']
namelist2= [' drhochemicalPotential',' chemicalPotential']
inputfile='../KortewegSources/kortewegCODEGEN.mpl'


subprocess.call( ['rm nskmaple.cc'],shell=True)
subprocess.call( ['maple '+inputfile],shell=True )



f=open( 'nskmaple.cc')
fnew=open( 'maplenew.cc','w')
flag = -1 
for line in f:
  if line[0]=='#':
    fnew.write('\n')
  else:
    newline1 = line.replace('theta','theta_')
    newline = newline1.replace('cst','cst_')
    for name in namelist:
      if newline.find( name ) != -1:
        flag=1
        break 
    for name2 in namelist2:
      if newline.find( name2 ) !=-1:
        flag=2
        break

    if flag==1:   
      newline = 'inline double'+name+' ( double rho ) const\n'
      fnew.write(newline)
      flag=-1
    elif flag==2:
      newline = 'inline double'+name2+' ( double rho ,double old ) const\n'
      fnew.write(newline)
      flag=-1
    else:
      fnew.write( newline )

subprocess.call( ['mv maplenew.cc ../KortewegSources/nskmaple.cc'], shell=True)

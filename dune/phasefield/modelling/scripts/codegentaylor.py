#! /usr/bin/env python
import os,subprocess, sys

folders = { 1:'../CoquelTaylorSources/' , 2:'../TaitTaylorSources/', 3:'../AltAltSources/'}

files = { 1:'coquelTaylorCODEGEN.mpl', 2:'taittaylorCODEGEN.mpl',3:'AltAltCODEGEN.mpl'}

namelist = [' helmholtz', ' pressure', ' a']
namelist2 = [ ' reactionSource',' dphireactionSource',' chemicalPotential',' dphichemicalPotential',' drhochemicalPotential']



number=int(sys.argv[1])

inputfile=folders[number]+files[number]

subprocess.call( ['rm maple.cc'],shell=True)
subprocess.call( ['maple '+inputfile],shell=True )

f=open( 'maple.cc')
fnew=open( 'maplenew.cc','w')
flag = -1 

for line in f:
  
  if line[0]=='#':
    
    fnew.write('\n')
  
  else:
    
    newline1 = line.replace('delta','delta_')
    newline2 = newline1.replace('beta', 'beta_') 
    newline3 = newline2.replace('A','A_')
    newline4 = newline3.replace('theta','theta_')
    newline = newline4.replace('alpha','alpha_')
    
    for name in namelist:
      if newline.find( name ) != -1:
       flag=1
       break 
    
    for name2 in namelist2:
      if newline.find( name2 ) != -1:
        flag=2
        break 
    
    if flag==1:   
      newline = 'inline double'+name+' ( double rho ,double phi ) const\n'
      fnew.write(newline)
      flag=-1
    elif flag==2:
      newline = 'inline double'+name2+' ( double rho ,double phi ,double old ) const\n'
      fnew.write(newline)
      flag=-1
    else:
      
      fnew.write( newline )

subprocess.call( ['mv maplenew.cc '+folders[number]+'/maple.cc'], shell=True)
subprocess.call( ['rm maple.cc'], shell=True )


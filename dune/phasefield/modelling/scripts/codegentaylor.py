#! /usr/bin/env python
import os,subprocess, sys

folders = { 1:'../CoquelTaylorSources/' , 2:'../CoquelTaylorSources/', 3:'../PhasefieldvanderWaalsSources/', 4:
'../CoquelTaylorSources/'}

files = { 1:'real' , 2:'coquelTaylor' , 3:'phasefieldvanderWaals', 4:'coquelTaylorM'}

namelist  = [' helmholtz', ' pressure', ' a']
namelist2 = [ ' reactionSource',' dphireactionSource',' drhoreactionSource',' chemicalPotential',' dphichemicalPotential',' drhochemicalPotential']
namelist3 = [' mwpliq', ' mwpvap', ' exactrho', ' exactv1', ' exactv2', ' exactphi',' exactsigma1',' exactsigma2', ' exactmu',' exacttau']
namelist4 = [' rhsRho', ' rhsV1', ' rhsV2', ' rhsPhi'] 


number=int(sys.argv[1])

inputfile=folders[number]+files[number]+'CODEGEN.mpl'

maplefile=files[number]+'maple.cc'
auxfile=files[number]+'mapleaux.cc'

subprocess.call( ['rm '+maplefile] , shell=True)
subprocess.call( ['maple '+inputfile] , shell=True )


rhofile=files[number]+'Rho.cc'
auxrhofile=files[number]+'Rhoaux.cc'

sourcefile=files[number]+'Source.cc'
auxsourcefile=files[number]+'Sourceaux.cc'



fSource=open( sourcefile )
fSourceNew=open( auxsourcefile, 'w')
flag=-1

for line in fSource:
    if line[0]=='#':
    
        fSourceNew.write('\n')
  
    else:
     
        l1= line.replace('delta','delta_')
        l2 = l1.replace('A','A_')
        l3 = l2.replace('visc1','mu1Liq_')
    
        for name in namelist4:
            if l3.find(name) !=-1:
                flag=1 
                break 
        
        if flag==1:
            newline='inline double'+name+' (double t, double x, double y ) const\n'
            fSourceNew.write( newline )
            flag=-1
        else:
            fSourceNew.write( l3 )

subprocess.call( ['mv '+auxsourcefile+' '+folders[number]+'/'+sourcefile], shell=True)
subprocess.call( ['rm '+sourcefile], shell=True )




frho=open( rhofile )
frhonew=open(  auxrhofile,'w' )
flag= -1


for line in frho:
    if line[0]=='#':
    
        frhonew.write('\n')
  
    else:
        l1= line.replace('delta','delta_')
        l2 = l1.replace('A','A_')
        l3 = l2.replace('visc1','mu1Liq_')
        
        for name in namelist3:
            if l3.find(name) !=-1:
                flag=1 
                break 
        
        if flag==1:
            newline='inline double'+name+' (double t, double x, double y ) const\n'
            frhonew.write( newline )
            flag=-1
        else:
            frhonew.write( l3 )


subprocess.call( ['mv '+auxrhofile+' '+folders[number]+'/'+rhofile], shell=True)
subprocess.call( ['rm '+rhofile], shell=True )

           
f=open( maplefile )
fnew=open( auxfile ,'w')
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

subprocess.call( ['mv '+auxfile+' '+folders[number]+'/'+maplefile], shell=True)
subprocess.call( ['rm '+maplefile], shell=True )



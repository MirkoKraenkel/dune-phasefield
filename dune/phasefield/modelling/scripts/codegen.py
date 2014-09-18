import os,subprocess, sys
folders = { 1:'../ConstRhoSources/' , 2:'../InvRhoSources/', 3:'../PhasefieldvanderWaalsSources/'}
files = { 1:'balancedCODEGEN.mpl', 2:'balancedhmodelCODEGEN.mpl', 3:'phasefieldvanderWaalsCODEGEN.mpl'}
namelist = [' helmholtz', ' reactionSource',' dphireactionSource',' chemicalPotential',' dphichemicalPotential',' drhochemicalPotential',' pressure', ' a']
number=int(sys.argv[1])
inputfile=folders[number]+files[number]
subprocess.call( ['rm maple.c'],shell=True)
subprocess.call( ['maple '+inputfile],shell=True )
f=open( 'maple.c')
fnew=open( 'maplenew.c','w')
flag = False
for line in f:
  if line[0]=='#':
    fnew.write('\n')
  else:
    newline1 = line.replace('delta','delta_')
    newline = newline1.replace('alpha','alpha_')
    for name in namelist:
      if newline.find( name ) != -1:
        flag=True
        break 
    if flag:   
      newline = 'inline double'+name+' ( double rho ,double phi ) const\n'
      fnew.write(newline)
      flag=False
    else:
      fnew.write( newline )
subprocess.call( ['mv maplenew.c '+folders[number]+'/maple.c'], shell=True)
subprocess.call( ['rm maple.c'], shell=True)

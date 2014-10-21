import os,subprocess, sys
folders = { 1:'../ConstRhoSources/' , 2:'../InvRhoSources/', 3:'../PhasefieldvanderWaalsSources/'}
files = { 1:'balancedCODEGEN.mpl', 2:'balancedhmodelCODEGEN.mpl', 3:'phasefieldvanderWaalsCODEGEN.mpl'}
namelist = [' helmholtz', ' reactionSource',' dphireactionSource',' chemicalPotential',' dphichemicalPotential',' drhochemicalPotential',' pressure', ' a']
namelist2 = [' rhosol', ' gradrho',' gradphi',' musol',' thetasol',' phiSource',' veloSource']
models={ 1:'balanced',2:'balancedh'}
number=int(sys.argv[1])
inputfile=folders[number]+files[number]
filename=models[number]+'.cc'
subprocess.call( ['rm '+filename],shell=True)
subprocess.call( ['rm maple.cc'],shell=True)
subprocess.call( ['maple '+inputfile],shell=True )
f=open( 'maple.cc')
fnew=open( 'maplenew.cc','w')
flag = False
for line in f:
  if line[0]=='#':
    fnew.write('\n')
  else:
    newline1 = line.replace('delta','delta_')
    newline2 = newline1.replace('beta', 'beta_') 
    newline3 = newline2.replace('A','A_')
    newline = newline3.replace('alpha','alpha_')
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
subprocess.call( ['mv maplenew.cc '+folders[number]+'/maple.cc'], shell=True)
subprocess.call( ['rm maple.cc'], shell=True )
if number==1 or number==2:
  f=open(filename)
  fnew=open( 'maplenew.cc','w')
  flag = False
  for line in f:
    if line[0]=='#':
      fnew.write('\n')
    else:
      newline1 = line.replace('delta','delta_')
      newline2 = newline1.replace('beta', 'beta_') 
      newline3 = newline2.replace('A','A_')
      newline = newline3.replace('alpha','alpha_')
      for name in namelist2:
        if newline.find( name ) != -1:
          flag=True
          break 
      if flag:   
        newline = 'inline double'+name+' ( double x ) const\n'
        fnew.write(newline)
        flag=False
      else:
        fnew.write( newline )
  fnew.close()
  f.close()
  subprocess.call( ['mv maplenew.cc ../../../../src/problems/mixedscheme/sourceprobCODEGEN/'+filename], shell=True)


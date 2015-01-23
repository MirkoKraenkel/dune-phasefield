import os,subprocess, sys
folders = { 1:'../Test/'}
files = { 1:'approxCODEGEN.mpl' }
namelist = [' helmholtz', ' reactionSource',' dphireactionSource',' chemicalPotential',' dphichemicalPotential',' drhochemicalPotential',' pressure', ' a']
namelist3= [' psi1',' psi2', ' mu1', ' mu2',' p1',' p2', ' dmu1', ' dmu2']
namelist2 = [' rhosol', ' gradrho',' gradphi',' musol',' thetasol',' phiSource',' veloSource']
models={ 1:'balanced',2:'balancedh',3:'vdW', 4:'freistuehler'}
number=int(sys.argv[1])
inputfile=folders[number]+files[number]
filename=models[number]+'.cc'
subprocess.call( ['rm '+filename],shell=True)
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
    newline5 = newline4.replace('ff1', 'psi1( rho )')
    newline6 = newline5.replace('ff2', 'psi2( rho )')
    newline7 = newline6.replace('m1', 'mu1( rho )')
    newline8 = newline7.replace('m0', 'mu2( rho )')
    newline9 = newline8.replace('ppp1', 'p1( rho )')
    newline10 = newline9.replace('ppp0', 'p2( rho )')
    newline11 = newline10.replace('dmu1', 'dmu1( rho )')
    newline12 = newline11.replace('dmu0', 'dmu2( rho )')
    newline = newline12.replace('alpha','alpha_')
    for name in namelist:
      if newline.find( name ) != -1:
       flag=1
       break 
    for nname in namelist3:
      if newline.find(nname) !=-1:
        flag=2
        break
    if flag==1:   
      print(name)
      newline = 'inline double'+name+' ( double rho ,double phi ) const\n'
      fnew.write(newline)
      flag=-1
    elif flag==2:
      newline = 'inline double'+nname+' ( double rho ) const\n'
      fnew.write(newline)
      flag=-1
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


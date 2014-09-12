import os,subprocess, sys
namelist = [' helmholtz',' chemicalPotential',' drhochemicalPotential',' pressure', ' a']
inputfile=sys.argv[1]
print(inputfile)
subprocess.call( ['rm maple.c'],shell=True)
subprocess.call( ['maple '+inputfile],shell=True )
f=open( 'maple.c')
fnew=open( 'maplenew.c','w')
flag = False
for line in f:
  if line[0]=='#':
    fnew.write('\n')
  else:
    newline1 = line.replace('theta','theta_')
    newline = newline1.replace('cst','cst_')
    for name in namelist:
      if newline.find( name ) != -1:
        flag=True
        break 
    if flag:   
      newline = 'inline double'+name+' ( double rho ) const\n'
      fnew.write(newline)
      flag=False
    else:
      fnew.write( newline )
subprocess.call( ['mv maplenew.c maple.c'], shell=True)

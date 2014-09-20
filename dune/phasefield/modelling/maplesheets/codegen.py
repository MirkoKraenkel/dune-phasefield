import os,subprocess, sys
namelist = [' rhosol', ' gradrho',' gradphi',' musol',' thetasol',' phiSource',' veloSource']

f=open( 'balancedh.c')
fnew=open( 'maplenew.c','w')
flag = False
for line in f:
  if line[0]=='#':
    fnew.write('\n')
  else:
    newline = line.replace('delta','delta_')
    for name in namelist:
      if newline.find( name ) != -1:
        flag=True
        break 
    if flag:   
      newline = 'inline double'+name+' ( double x) const\n'
      fnew.write(newline)
      flag=False
    else:
      fnew.write( newline )
fnew.close()
f.close()
subprocess.call( ['mv maplenew.c balancedh.c'], shell=True)

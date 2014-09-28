import os,subprocess, sys
namelist = [' rhosol', ' gradrho',' gradphi',' musol',' thetasol2',' thetasol2',' phiSource',' veloSource']
models=['balanced','balancedh']
number=int(sys.argv[1])
filename=models[number]+'.c'
f=open(filename)
fnew=open( 'maplenew.c','w')
flag = False
for line in f:
  if line[0]=='#':
    fnew.write('\n')
  else:
    newline2 = line.replace('delta','delta_')
    newline=newline2.replace('lambda','lambda_')
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
subprocess.call( ['mv maplenew.c '+filename], shell=True)

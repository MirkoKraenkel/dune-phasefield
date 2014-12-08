#! /usr/bin/env python
import subprocess,sys,time
order=sys.argv[1]
order='_P'+order
if len(sys.argv)==3:
  tag=sys.argv[2]
  tag='_'+tag
else:
  tag=''
myday=time.strftime("%d_%m_%Y")
print(myday)
programms=['phasefield_coupling']
for prog in programms:
  name=prog+order
  callstring='./'+name+' parameter path:data'+tag+order+'_'+myday+'/'+name
  print(callstring)
  subprocess.call([callstring],shell=True)
  

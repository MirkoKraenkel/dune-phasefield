import os, subprocess, sys
rho='rho'
v='v'
vx='vx'
vy='vy'
u='u'
phi='phi'
mu='mu'
tau='tau'
sigmax='sigmax'
sigmay='sigmay'
###################################################################################################################

acvariables1d={ phi : 0 , tau : 1 , sigmax : 2 }
acvarlist=list(acvariables1d.keys())
acvarlist.sort()


acelemcouplings1d=[ 
            ( phi,phi),(phi,tau),
            (tau,phi),(tau,tau),(tau,sigmax),
            (sigmax,phi),(sigmax,sigmax)
]
acinterseccouplings1d=[  
            ( phi,phi),
            (tau,phi),(tau,sigmax),
            (sigmax,phi)
]
f=open('acelCouplings1d_CODEGEN.c','w')
i=0
for cp in acelemcouplings1d:
  v1=acvariables1d[ cp[0] ]
  v2=acvariables1d[ cp[1] ]
  newline='elementCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
f=open('acisCouplings1d_CODEGEN.c','w')
i=0
for cp in acinterseccouplings1d:
  v1=acvariables1d[ cp[0] ]
  v2=acvariables1d[ cp[1] ]
  newline='intersectionCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
###################################################################################################################################
nvstvariables1d={ rho : 0 , v : 1 , mu : 2  }
nvstvarlist=list(nvstvariables1d.keys())
nvstvarlist.sort()


nvstelemcouplings1d=[ ( rho , rho ), ( rho , v ), 
            ( v , rho ) ,( v , v), ( v , mu ),
            ( mu, rho),(mu,v),(mu,mu),
]
nvstinterseccouplings1d=[  ( rho , rho ), ( rho , v ),( rho , mu),
            ( v , rho ) ,( v , v),( v , mu )
]
f=open('nvstelCouplings1d_CODEGEN.c','w')
i=0
for cp in nvstelemcouplings1d:
  v1=nvstvariables1d[ cp[0] ]
  v2=nvstvariables1d[ cp[1] ]
  newline='elementCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
f=open('nvstisCouplings1d_CODEGEN.c','w')
i=0
for cp in nvstinterseccouplings1d:
  v1=nvstvariables1d[ cp[0] ]
  v2=nvstvariables1d[ cp[1] ]
  newline='intersectionCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
#################################################################################################################################

acvariables2d={ phi : 0 ,  tau : 1 , sigmax : 2 , sigmay : 3 }
acvarlist=list(acvariables2d.keys())
acvarlist.sort()

acelemcouplings2d=[( phi,phi),(phi,tau),
            (tau,phi),(tau,tau),(tau,sigmax),(tau,sigmay),
            (sigmax,phi),(sigmax,sigmax),
            (sigmay,phi),(sigmay,sigmay)
]
acinterseccouplings2d=[ ( phi,phi),
            (tau,phi),(tau,sigmax),(tau,sigmay),
            (sigmax,phi),
            (sigmay,phi)
]

f=open('acelCouplings2d_CODEGEN.c','w')
i=0
for cp in acelemcouplings2d:
  v1=acvariables2d[ cp[0] ]
  v2=acvariables2d[ cp[1] ]
  newline='elementCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
f=open('acisCouplings2d_CODEGEN.c','w')
i=0
for cp in acinterseccouplings2d:
  v1=acvariables2d[ cp[0] ]
  v2=acvariables2d[ cp[1] ]
  newline='intersectionCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()

######################################################################################################################
nvstvariables2d={ rho : 0 , vx : 1 , vy : 2 , mu : 3  }
nvstvarlist=list(nvstvariables2d.keys())
nvstvarlist.sort()

nvstelemcouplings2d=[ ( rho , rho ),( rho , vx ),(rho,vy),
             (vx , rho ),( vx , vx ),( vx , vy ), ( vx , mu ),
            ( vy , rho ),( vy , vx ),( vy , vy ), ( vy , mu ),
            ( mu, rho),(mu,vx),(mu,vy),(mu,mu)
]
nvstinterseccouplings2d=[  ( rho , rho ), ( rho , vx ),( rho , vy ),( rho , mu),
            ( vx , rho ),( vx , vx ),(vx,vy), ( vx , mu ),
            ( vy , rho ),( vy , vx ),(vy,vy), ( vy , mu )
]

f=open('nvstelCouplings2d_CODEGEN.c','w')
i=0
for cp in nvstelemcouplings2d:
  v1=nvstvariables2d[ cp[0] ]
  v2=nvstvariables2d[ cp[1] ]
  newline='elementCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
f=open('nvstisCouplings2d_CODEGEN.c','w')
i=0
for cp in nvstinterseccouplings2d:
  v1=nvstvariables2d[ cp[0] ]
  v2=nvstvariables2d[ cp[1] ]
  newline='intersectionCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
########################################################################################################################

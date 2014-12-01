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


variables1d={ rho: 0 , v: 1 ,phi: 2, mu : 3 , tau : 4, sigmax:  5 }
varlist=list(variables1d.keys())
varlist.sort()

elemcouplings1d=[ ( rho , rho ), ( rho , v ), 
            ( v , rho ) ,( v , v), ( v, phi ),( v , mu ),( v , tau ),
            ( phi , rho ),( phi , v ), ( phi,phi),(phi,tau),
            ( mu, rho),(mu,v),(mu,phi),(mu,mu),
            ( tau,rho),(tau,phi),(tau,tau),(tau,sigmax),
            (sigmax,phi),(sigmax,sigmax)
]
interseccouplings1d=[  ( rho , rho ), ( rho , v ),( rho , mu),
            ( v , rho ) ,( v , v), ( v, phi ),( v , mu ),( v , tau ),
           ( phi , v ), ( phi,phi),
            (tau,phi),(tau,sigmax),
            (sigmax,phi)
]
f=open('elementCouplings1d_CODEGEN.c','w')
i=0
for cp in elemcouplings1d:
  v1=variables1d[ cp[0] ]
  v2=variables1d[ cp[1] ]
  newline='elementCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
f=open('intersectionCouplings1d_CODEGEN.c','w')
i=0
for cp in interseccouplings1d:
  v1=variables1d[ cp[0] ]
  v2=variables1d[ cp[1] ]
  newline='intersectionCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()

variables2d={ rho: 0 , vx: 1 , vy:2 ,phi: 3, mu : 4 , tau : 5, sigmax:  6 ,sigmay: 7}
varlist=list(variables2d.keys())
varlist.sort()

elemcouplings2d=[ ( rho , rho ),( rho , vx ),(rho,vy),
             (vx , rho ),( vx , vx ),( vx , vy ), ( vx, phi ),( vx , mu ),( vx , tau ),
            ( vy , rho ),( vy , vx ),( vy , vy ), ( vy, phi ),( vy , mu ),( vy , tau ),
            ( phi , rho ),( phi , vx ),( phi , vy ), ( phi,phi),(phi,tau),
            ( mu, rho),(mu,vx),(mu,vy),(mu,phi),(mu,mu),
            ( tau,rho),(tau,phi),(tau,tau),(tau,sigmax),(tau,sigmay),
            (sigmax,phi),(sigmax,sigmax),
            (sigmay,phi),(sigmay,sigmay)
]
interseccouplings2d=[  ( rho , rho ), ( rho , vx ),( rho , vy ),( rho , mu),
            ( vx , rho ),( vx , vx ),(vx,vy), ( vx, phi ),( vx , mu ),( vx , tau ),
            ( vy , rho ),( vy , vx ),(vy,vy), ( vy, phi ),( vy , mu ),( vy , tau ),
            ( phi , vx ),( phi , vy ), ( phi,phi),
            (tau,phi),(tau,sigmax),(tau,sigmay),
            (sigmax,phi),
            (sigmay,phi)
]
f=open('elementCouplings2d_CODEGEN.c','w')
i=0
for cp in elemcouplings2d:
  v1=variables2d[ cp[0] ]
  v2=variables2d[ cp[1] ]
  newline='elementCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
f=open('intersectionCouplings2d_CODEGEN.c','w')
i=0
for cp in interseccouplings2d:
  v1=variables2d[ cp[0] ]
  v2=variables2d[ cp[1] ]
  newline='intersectionCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()


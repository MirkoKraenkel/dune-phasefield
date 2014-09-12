import os, subprocess, sys
rho='rho'
v='v'
u='u'
phi='phi'
mu='mu'
tau='tau'
sigmax='sigmax'
sigmay='sigmay'


f=open('elementCouplings1d_CODEGEN.c','w')
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
interseccouplings1d=[  ( rho , rho ), ( rho , v ), 
            ( v , rho ) ,( v , v), ( v, phi ),( v , mu ),( v , tau ),
           ( phi , v ), ( phi,phi),
            (tau,tau),(tau,sigmax),<F2>
            (sigmax,phi),(sigma,sigmax)
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


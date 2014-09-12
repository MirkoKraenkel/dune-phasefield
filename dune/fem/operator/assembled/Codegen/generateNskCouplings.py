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
variables1d={ rho: 0 , v: 1 ,mu : 2 , sigmax: 3 }
varlist=list(variables1d.keys())
varlist.sort()
elemcouplings1d=[ ( rho , rho ), ( rho , v ), (rho,mu),
            ( v , rho ) ,( v , v), ( v , mu ),
            ( mu, rho),(mu,v),(mu,mu),(mu,sigmax),
            (sigmax,rho),(sigmax,sigmax)
]
interseccouplings1d=[  ( rho , rho ), ( rho , v ), (rho, mu),
            ( v , rho ) ,( v , v), ( v , mu ),
            (mu ,mu ),( mu ,sigmax),
            (sigmax,rho),(sigmax,sigmax)
]
variables2d={ rho: 0 , v: 1 ,u : 2, mu:3 , sigmax: 4, sigmay:5 }
varlist2d=list(variables2d.keys())
varlist2d.sort()
elemcouplings2d=[ ( rho , rho ) , ( rho , v ), (rho,u), (rho,mu),
                  ( v , rho ) , ( v , v ) , ( v , u ) , ( v , mu ),
                  ( u , rho ) , ( u , v ) , ( u , u ) , ( u , mu ),
                  ( mu, rho ) , ( mu ,v ) , ( mu , u ) , ( mu , mu ) , ( mu , sigmax ) , (  mu  , sigmay ),
                  ( sigmax , rho ) , ( sigmax , sigmax ),
                  ( sigmay , rho ) , ( sigmay , sigmay ) ]
interseccouplings2d=[ ( rho , rho ) , ( rho , v ), (rho,u), (rho,mu),
                      ( v , rho ) , ( v , v ) , ( v , u ) , ( v , mu ),
                      ( u , rho ) , ( u , v ) , ( u , u ) , ( u , mu ),
                      ( mu , mu ) , ( mu , sigmax ) , (  mu  , sigmay ),
                      ( sigmax , rho ) , ( sigmax , sigmax ),
                      ( sigmay , rho ) , ( sigmay , sigmay ) ]

f=open('nsk_element_1d_CODEGEN.c','w')

i=0
for cp in elemcouplings1d:
  v1=variables1d[ cp[0] ]
  v2=variables1d[ cp[1] ]
  newline='elementCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
f=open('nsk_intersection_1d_CODEGEN.c','w')
i=0
for cp in interseccouplings1d:
  v1=variables1d[ cp[0] ]
  v2=variables1d[ cp[1] ]
  newline='intersectionCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
f=open('nsk_element_2d_CODEGEN.c','w')
i=0
for cp in elemcouplings2d:
  v1=variables2d[ cp[0] ]
  v2=variables2d[ cp[1] ]
  newline='elementCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()
f=open('nsk_intersection_2d_CODEGEN.c','w')
i=0
for cp in interseccouplings2d:
  v1=variables2d[ cp[0] ]
  v2=variables2d[ cp[1] ]
  newline='intersectionCouplings_[ '+str( i )+' ] = std::make_pair('+str(v1)+','+str(v2)+');//('+cp[0]+','+cp[1]+')\n'
  f.write(newline)
  i=i+1
f.close()


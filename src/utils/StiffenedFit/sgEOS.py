#! /usr/bin/env python
import sys,pickle

Tsat=517.15
    
Psat=8.5879
    
rhoSat=46.166

energy=2563.6

enthalp=2749.6

entrop=5.79059

cP=6.2187

sndSpd=480.73
gamma=sndSpd*sndSpd/(cP*Tsat*1000)+1
print(gamma)
estar=enthalp-cP*Tsat
print(estar)
pInf=1/gamma*( (gamma-1)*rhoSat*(energy-estar)*1000-Psat*1000000 )
print( pInf )
p2=((gamma-1)*rhoSat*(energy-estar)*1000-gamma*pInf)*1e-6
print(p2)
cV=(Psat*1e+6+pInf)/(rhoSat*Tsat*(gamma-1))
print(cV)

a=cV*Tsat*(gamma-1)
print(a*rhoSat-pInf)*1e-6

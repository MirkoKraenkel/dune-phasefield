#! /usr/bin/env python
import sys,pickle
from pylab import *
import math
import numericalunits  as nu
import numpy as np
#import matplotlib as plt
nu.reset_units()

Tcrit   = 647.096*nu.K

Pcrit   = 22.064*nu.MPa


rhoCrit =322*3*nu.kg/(nu.m*nu.m*nu.m)

lam=nu.J/nu.m

L = 1e-4*nu.m
g = 9.81*nu.m/(nu.s*nu.s)

velgrav=sqrt(L*g)
print("gravitationalVelocity= "+str(velgrav*nu.s/nu.m))
print('gravitationalTime = '+str(L/velgrav/nu.s))
eref     = Pcrit/rhoCrit
visc  = 9.4e-5*nu.Pa*nu.s

criticalTime=sqrt(L*L*rhoCrit/Pcrit)
print('criticalTime = '+str(criticalTime/nu.s))

criticalVelocity=sqrt(Pcrit/rhoCrit)
print('criticalVelocity= '+str(criticalVelocity*nu.s/nu.m))


print('Eref='+str( eref*rhoCrit /( (nu.J/nu.g)*(nu.kg/(nu.m*nu.m*nu.m)))))
# [ Tsat , Psat , rhoSat , energy ,enthalp , entrop , cP , sndSpd  ]
#valVapour=[ 573.15 , 8.5879 , 46.168, 2563.6 , 2749.6 , 5.7059 , 6.2197 , 480.73]
#valWater=[ 573.15 , 8.5879 , 712.14 , 1332.9 , 1345.0 , 3.2552 , 5.7504 , 909.4 ]
valVapour=[ 550.316,6.12013, 31.4907,2588.92,2783.27,5.88070,4.93455,493.296]
valWater=[ 550.316,6.12013,755.75,1212.54,1220.64,3.03972,5.23362,   1027.70]

deltaRho=(valWater[2]-valVapour[2])*nu.kg/(nu.m*nu.m*nu.m)
#deltaRho/=rhoCrit
#print("deltaRho= "+str(deltaRho))

Cahn=1e-3

#vals=[valVapour,valWater]
##print(valVapour[0]*nu.K/Tcrit)
erefbulk = (eref-valVapour[1])/rhoCrit

sigma    = 0.0197225*nu.J/(nu.m*nu.m)

#Weber=rhoCrit*velgrav*velgrav*L/sigma
Weber=Pcrit*L/sigma
Weber2=rhoCrit*criticalVelocity*criticalVelocity*L/sigma

Froud=(criticalTime/criticalVelocity)*deltaRho/rhoCrit*g

Eotvos=(valWater[2]-valVapour[2])*(nu.kg/(nu.m*nu.m*nu.m))*g*L*L/sigma
print("Eoetvoes Number= "+str(Eotvos)+'\n')

print(' Critical Weber Number = '+str(Weber2)+'; 1/We = '+str(1/Weber2)+'\n')
print(' Rhs = '+str( Froud )+'\n')

print(' WFactor= '+str(1/(Cahn*Weber2))+' GradFactor = '+str(Cahn/Weber2)+'\n')

#Reynolds=rhoCrit*velgrav*L/visc
Reynolds=sqrt(Pcrit*rhoCrit)*L/visc
Re2=rhoCrit*criticalVelocity*L/visc

print(' Critical Reynolds Number = '+str(Re2)+'\n')
print(' Critical Capillary = '+str(Weber2/Re2)+'\n')

#print(' gamma NSK(=1/We^2) = '+str(1/(Weber2*Weber2))+'\n')

print(' Gravity Reynolds Number'+str(Reynolds)+'\n')
#print('Gravity Capillary = '+str(Weber/Reynolds)+'\n')

class StiffenedEnergy:
    
    def __init__( self , vals ):
        self.Tsat    = vals[0]*nu.K
        self.Tsat   /= Tcrit
        self.Psat    = vals[1]*nu.MPa
        self.Psat   /= Pcrit
        self.rhoSat  = vals[2]*nu.kg/(nu.m*nu.m*nu.m)
        self.rhoSat /= rhoCrit
        self.energy  = vals[3]*nu.kJ/nu.kg
        self.energy /= eref
        self.enthalp = vals[4]*nu.kJ/nu.kg
        self.enthalp/= eref
        self.entrop  = vals[5]*nu.J/(nu.g*nu.K)
        self.entrop /= eref/Tcrit
        self.cP      = vals[6]*nu.J/(nu.g*nu.K)
        self.cP     /= eref/Tcrit
        self.sndSpd  = vals[7]*(nu.m/nu.s)
        self.sndSpd /= criticalVelocity 
        self.gamma=self.sndSpd*self.sndSpd/(self.cP*self.Tsat)+1
        self.estar=self.enthalp-self.cP*self.Tsat
        self.pInf=1/self.gamma*( (self.gamma-1)*self.rhoSat*(self.energy-self.estar)-self.Psat )
        self.cV=(self.Psat+self.pInf)/(self.rhoSat*self.Tsat*(self.gamma-1))
        self.cpot=self.enthalp-self.Tsat*self.entrop
        self.a=self.cV*self.Tsat*(self.gamma-1)
        
        self.b=self.cV*self.Tsat-self.cV*self.Tsat*(self.gamma-1)*np.log(self.rhoSat)+self.cV*self.Tsat*(self.gamma-1)-self.entrop*self.Tsat+self.estar
        print(" soundspeed= "+str(self.sndSpd)+"\n") 
    def Afree(self,rho):
    
        Afree=self.cV*self.Tsat*(1+(self.gamma-1)*np.log(rho/self.rhoSat))-self.entrop*self.Tsat+self.pInf/rho+self.estar
        
        Afree*=rho
        
        return Afree 
    
    def A2(self,rho):
        helm=self.a*rho*np.log(rho)+(self.b-self.a)*rho+self.pInf
        return helm
    def pressure(self,rho):
        pres=(self.gamma-1)*self.cV*self.Tsat*rho-self.pInf 
        return pres   
    def pressure2(self,rho):
        pres=self.a*rho-self.pInf
        return pres
    
    def mu(self,rho):
        m=self.a*np.log(rho)+self.b*rho
        return m

if __name__ == '__main__':
    vapA = StiffenedEnergy( valVapour )
    liqA = StiffenedEnergy( valWater )
    
#    print('a_Vap= '+ str(vapA.a ) )
#    print('b_Vap= '+ str(vapA.b ) )
#    print('p_Vap= '+ str(vapA.pInf) )
#    print('rhoVap= '+str(vapA.rhoSat) )
#    print('a_Liq= '+ str(liqA.a ) )
#    print('b_Lip= '+ str(liqA.b #) #)
#    print('p_Liqp= '+ str(liqA.pInf) )
#    print('rhoLiq= '+str(liqA.rhoSat) )
#    print('Vappress= '+str(vapA.pressure2( vapA.rhoSat) ))
#    print('Liqpress= '+str(liqA.pressure2( liqA.rhoSat) ))
#    print('Vapmu= '+str(vapA.mu( vapA.rhoSat) ))
#    print('Liqmu= '+str(liqA.mu( liqA.rhoSat) ))
#    print( str(vapA.enthalp)+' - '+str(vapA.Tsat)+' * '+str(vapA.entrop))
#    print('Vapcpot= '+str(vapA.cpot) )
#    print( str(liqA.enthalp)+' - '+str(liqA.Tsat)+' * '+str(liqA.entrop))
#    print('Liqcpot= '+str(liqA.cpot) )


    rr = np.arange(0.001, 1,1e-3)
    fvap = vapA.A2( rr )
    fliq = liqA.A2( rr )


    pliq = liqA.pressure( rr )
    
    filename='paramWater_T_'+str(vapA.Tsat)
    f=open(filename,'w')
    
    f.write('sg.bl: '+str( liqA.b)+'\n')
    f.write('sg.bv: '+str( vapA.b)+'\n')
    f.write('sg.cl: '+str( liqA.a)+'\n')
    f.write('sg.cv: '+str( vapA.a)+'\n')
    f.write('sg.dl: '+str( liqA.pInf)+'\n')
    f.write('sg.dv: '+str( vapA.pInf)+'\n')
    f.write('phasefield.mwpvap: '+str( vapA.rhoSat )+'\n')
    f.write('phasefield.mwpliq: '+str( liqA.rhoSat )+'\n')
    f.close()

#    plot(rr,fvap,rr, fliq)
#    grid(True)
#    show()

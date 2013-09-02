#!/bin/sh
#make clean
make POLORDER=3 phasefield_nc ; mv phasefield_nc nc3
make clean
make POLORDER=2 phasefield_nc ; mv phasefield_nc nc2
make clean
make POLORDER=1 phasefield_nc ; mv phasefield_nc nc1
#make clean
#make POLORDER=0 phasefield_nc ; mv phasefield_nc nc0




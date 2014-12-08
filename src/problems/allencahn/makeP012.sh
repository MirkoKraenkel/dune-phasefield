#!/bin/sh
#make clean
#make POLORDER=3 phasefield_nc ; mv phasefield_nc nc3
make clean
make POLORDER=2 allencahn ; mv allencahn ac2
make clean
make POLORDER=1 allencahn ; mv allencahn ac1
make clean
make POLORDER=0 allencahn ; mv allencahn ac0




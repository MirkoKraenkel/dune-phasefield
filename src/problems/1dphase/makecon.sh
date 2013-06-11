#!/bin/sh
make clean
make POLORDER=2 phasefield_nc ; mv phasefield_nc nc2
make clean
make POLORDER=1 phasefield_nc ; mv phasefield_nc nc1
make clean
make POLORDER=0 phasefield_nc ; mv phasefield_nc nc0


make clean
make POLORDER=2 phasefield ; mv phasefield con2
make clean
make POLORDER=1 phasefield ; mv phasefield con1
make clean
make POLORDER=1 phasefield ; mv phasefield con0


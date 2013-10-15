#!/bin/sh
#make clean
#make GRIDDIM=2 phasefield_nc ; mv phasefield_nc nc_2d_nothetavisc
make clean
make GRIDDIM=1 phasefield_nc ; mv phasefield_nc nc_old
make clean
make GRIDDIM=1 phasefield_nctr mv phasefield_nctr nc_new



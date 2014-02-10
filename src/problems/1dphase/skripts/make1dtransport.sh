#!/bin/sh
#make clean
#make GRIDDIM=2 phasefield_nc ; mv phasefield_nc nc_2d_nothetavisc
rm trans_old
rm trans_new
make clean
make GRIDTYPE=SPGRID GRIDDIM=1 phasefield_nc ; mv phasefield_nc trans_old
make clean
make GRIDTYPE=SPGRID GRIDDIM=1 phasefield_nctr ;mv phasefield_nctr trans_new



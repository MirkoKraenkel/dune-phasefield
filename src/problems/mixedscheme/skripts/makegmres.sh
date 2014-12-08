#!/bin/sh
x=2
MAXORD=$1
echo $MAXORD
while [ $x -le $MAXORD ]
do
make clean
make POLORDER=$x phasefielddggmres ; mv phasefielddggmres pfP$x$HOSTNAME$2
x=$(( $x + 1))
done

#make clean
#make POLORDER=1 phasefield_nc ; mv phasefield_nc nc1
#make clean
#make POLORDER=0 phasefield_nc ; mv phasefield_nc nc0


#make clean
#make POLORDER=2 phasefield ; mv phasefield con2
#make clean
#make POLORDER=1 phasefield ; mv phasefield con1
#make clean
#make POLORDER=1 phasefield ; mv phasefield con0


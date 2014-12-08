#!/bin/sh
x=1
MAXORD=$1
d=0.01
l=4
echo $MAXORD
while [ $x -le $MAXORD ]
do
echo $d
./phasefielddggmres scheme:2d$x  phasefield.delta:$d fem.adaptation.finestLevel:$((2*$x+$l)) 
d=`echo "scale=10;$d/4 "  | bc`
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


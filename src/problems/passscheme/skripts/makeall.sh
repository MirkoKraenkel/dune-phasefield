#!/bin/sh
x=0
if [$# -lt 2]
then
  echo "Usage: PolOrder Schmeme: (1=conservative, 2=nocncon1, 3=noncton1)"
  exit
fi
MAXORD=$1
echo $MAXORD

case $2 in
1) SCHEMENAME=phasefield
    ;;
2) SCHEMENAME=phasefied_nc
    ;;
3) SCHEMENAME=phasefield_nctr
    ;;
*) echo "Wrong Scheme number!"
  
while [ $x -le $MAXORD ]
do
make clean
make POLORDER=$x $(SCHEMENAME) ; mv $(SCHEMENAME) $(SCHEMENAME)_P$x
x=$(( $x + 1))
done



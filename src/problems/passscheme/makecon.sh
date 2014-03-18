#!/bin/sh

make clean
make POLORDER=2 phasefield ; mv phasefield con2
make clean
make POLORDER=1 phasefield ; mv phasefield con1
make clean
make POLORDER=1 phasefield ; mv phasefield con0


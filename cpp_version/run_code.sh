#!/bin/bash

gcan=1.0
gKv1=1.0
gKCa=0.25
Imin=0
Imax=20
Ib=5

g++ -std=c++20 -O2 -I ../boost main.cpp -o neuron 
./neuron -I $Imin $Ib $Imax -can $gcan -ka $gKv1 -kca $gKCa 
rm neuron
mv dat0 current_ramp_up_data.txt
mv dat1 current_tent_data.txt
rm dat2
mv dat3 current_step_data.txt
julia plot.jl
open plot.png
rm *.txt
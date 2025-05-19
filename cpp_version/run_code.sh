#!/bin/bash

gcan=1.0
gKv1=1.0
gKCa=0.25
Imin=0
Imax=20
Ib=5

# Using OMP on Mac can be annoying. This setup uses the homebrew g++ compiler and the homebrew install of OMP.
# The main.cpp file is set up to conditionally import OMP, so you can also just remove the -fopenmp flags when compiling
# and everything work normally (but not using OMP)
if [[ "$OSTYPE" == "darwin"* ]]; then
    g++-14 -std=c++20 -O2 -fopenmp -I ../boost/ main.cpp -o neuron 
else
    g++ -std=c++20 -O2 -fopenmp -I ../boost main.cpp -o neuron 
fi
./neuron -I $Imin $Ib $Imax -can $gcan -ka $gKv1 -kca $gKCa 
rm neuron
julia plot.jl
open plot.png
rm *.txt
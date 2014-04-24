#!/bin/bash

for no in 4 8 16 32
do
	for N in 4 8 16 32
	do
		matlab -nosplash -nodesktop -nodisplay -r "exp_lorenz_enkf(128,$N,$no);quit" >log/lorenz_$N_$no.log &
	done
done	
#!/bin/bash

ts=1
for no in 8
do
	for N in 4 8 16
	do
		matlab -nosplash -nodesktop -nodisplay -r "exp_swe_enkf($N,,$no,$ts);quit" >"log/lorenz_$N_$no_$ts.log" &
	done
done	
#!/bin/bash

ts=5
for no in 4 8 16 32
do
	for N in 4 8 16
	do
		matlab -nosplash -nodesktop -nodisplay -r "exp_swe_enkf($N,$no,$ts);quit" >"log/swe_$N_$no_$ts.log" &
	done
done	
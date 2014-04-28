#!/bin/bash

ts=5
for no in 4 8 16 
do
	for N in 4 8 16
	do
		matlab -nosplash -nodesktop -nodisplay -r "exp_swe_enkf($N,$ts,$no);quit" >"log/swe_${N}_${no}_${ts}.log" &
	done
done	
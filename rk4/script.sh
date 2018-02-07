#!/bin/bash

# params
LOGFILE=log.txt
N_RANGE="101 501 751"

for N in ${N_RANGE}
do
	echo ${N} 
	./bin/repr_rk45 -n ${N} >> ${LOGFILE}
	rm -rf ./bin/*.csv
done



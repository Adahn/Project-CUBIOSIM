#!/bin/bash

# params
LOGFILE=log.txt
N_RANGE="101 501 1001 5001 10001 15001 20001 50001 100001"

for N in ${N_RANGE}
do
	echo ${N} 
	./bin/repr_rk4 -n ${N} >> ${LOGFILE}
	rm -rf ./bin/*.csv
done



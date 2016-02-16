#!/bin/bash

ARCH=x86
INIT=8388607
RANGE=16777216
REP=5

mkdir -p results

for file in results/*.csv; do mv "$file" "${file}.old"; done
for file in results/*.dat; do mv "$file" "${file}.old"; done

for prog in "BSTTK" "GreenBST" "LFBST"  "SVEB"  "citrus" "CBTree"  
do
	for upd in 0 50
	do
		for thread in 1 6 12 18 24
		do

			DELTA=()
			MEAN=()
			M2=()
			LABEL=()
			STDEV=()

			#echo "#D: $prog - $thread - $upd" >> results/$prog.$ARCH.csv

			n=1

			STR="5000000,${thread},${upd}"

			while [ "$n" -le "$REP" ]
			do
				../$prog/$prog.energy -s 0 -n $thread -u $upd -i $INIT -r $RANGE > $prog.$ARCH 2>&1
				while [ $? -eq 99 ]; do
					../$prog/$prog.energy -s 0 -n $thread -u $upd -i $INIT -r $RANGE > $prog.$ARCH 2>&1
				done
				echo $(grep "0:" $prog.$ARCH) >> results/$prog-$INIT.$n.dat
				echo $(grep "#P" $prog.$ARCH) >> results/$prog-$INIT.$n.dat
				
				STATS=$(grep "#D" $prog.$ARCH)
				
				COUNTERS=0
				while read -r line; do
					#echo $line

					if [ "${LABEL[$COUNTERS]}" = "" ]; then
						LABEL[$COUNTERS]=$(echo $line | cut -d',' -f2)
						DELTA[$COUNTERS]=0
						MEAN[$COUNTERS]=0
						M2[$COUNTERS]=0	
					fi

					X=$(echo $line | cut -d',' -f3)
					#echo "X=$X"
					
					#echo "$X - ${MEAN[$COUNTERS]}"
					DELTA[$COUNTERS]=$(echo "scale=10; $X - ${MEAN[$COUNTERS]}" | bc)
					#echo "DELTA=${DELTA[$COUNTERS]}"
						
					#echo "${MEAN[$COUNTERS]} + (${DELTA[$COUNTERS]}/$n)"
					MEAN[$COUNTERS]=$(echo "scale=10; ${MEAN[$COUNTERS]} + (${DELTA[$COUNTERS]}/$n)" | bc)
					#echo "MEAN=${MEAN[$COUNTERS]}"

					#echo "${M2[$COUNTERS]} + (${DELTA[$COUNTERS]}*($X - ${MEAN[$COUNTERS]}))"
					M2[$COUNTERS]=$(echo "scale=10; ${M2[$COUNTERS]} + (${DELTA[$COUNTERS]}*($X - ${MEAN[$COUNTERS]}))" | bc)
					#echo "M2=${M2[$COUNTERS]}"

					if [ $n == $REP ]; then
						STDEV[$COUNTERS]=$(echo "scale=10; sqrt(${M2[$COUNTERS]}/($n-1))" | bc)
						#echo ${LABEL[$COUNTERS]},${MEAN[$COUNTERS]},${STDEV[$COUNTERS]} >> results/$prog.$ARCH.csv
						STR="${STR},${MEAN[$COUNTERS]},${STDEV[$COUNTERS]}"
					fi

					COUNTERS=$((COUNTERS+1))

				done <<< "$STATS"
				n=$((n+1))
			done

			echo $STR >> results/$prog.$ARCH.csv

		done

	done
done


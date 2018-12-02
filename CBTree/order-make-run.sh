#! /bin/bash

init=16777215
rang=33554432

thrd=36

jemalloc=`/opt/jemalloc/bin/jemalloc-config --libdir`/libjemalloc.so.`/opt/jemalloc/bin/jemalloc-config --revision`
bin=CBTree

st=5
ed=12
stamp=`date '+%y%m%d%H%M%S'`
for i in $(eval echo {$st..$ed})
do 
	ORDER=$((2**$i)) 
	ORDER=${ORDER} make withorder ARCH=x86_64

	for UPD in 10 50
	do
		./${bin}.energy.${ORDER} -n $thrd -u $UPD -i $init -r $rang >> $bin.order.stdmalloc.$stamp 2>&1
		LD_PRELOAD=$jemalloc ./${bin}.energy.${ORDER} -n $thrd -u $UPD -i $init -r $rang >> $bin.order.jemalloc.$stamp 2>&1
	done
done

echo "done. Time: $stamp"


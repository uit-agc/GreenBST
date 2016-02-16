#!/bin/bash

for prog in "BSTTK" "GreenBST" "LFBST"  "SVEB"  "citrus" "CBTree"  
do
	cd ../$prog
	make clean
	make
done

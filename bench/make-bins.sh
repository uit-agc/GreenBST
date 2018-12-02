#!/bin/bash

for prog in "BSTTK" "GreenBST" "LFBST" "SVEB" "citrus" "CBTree" "abtree" "bwtree"
do
	cd ../$prog
	make clean
	make
done


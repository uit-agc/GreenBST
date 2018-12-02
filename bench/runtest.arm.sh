#!/bin/bash

ARCH=ARM
SUFFIX=energy

INIT=1118481
RANGE=2236964

UPD=(10 50)
THREAD=(1 2 3 4)
PROGS=("BSTTK" "GreenBST" "LFBST" "SVEB" "citrus" "CBTree" "abtree" "bwtree")

REP=5

# Run bench
my_dir="$(dirname "$0")"
. "$my_dir/runtest.in.sh"

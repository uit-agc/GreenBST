#!/bin/bash

ARCH=x86
SUFFIX=energy

INIT=16777215
RANGE=33554432

UPD=(10 50)
THREAD=(1 9 18 27 36)
PROGS=("BSTTK" "GreenBST" "LFBST" "SVEB" "citrus" "CBTree" "abtree" "bwtree")

REP=5

# Run bench
my_dir="$(dirname "$0")"
. "$my_dir/runtest.in.sh"

#!/bin/bash

# This script compiles and runs the simulation
# How to use it: 
# run.sh -d <simulation-directory> -n <number-of-cores> -c <compilation-type> 

BIN="simulate"
DIRPREFIX="simulation-"
AUX=`ls simulations | grep $DIRPREFIX | tail -1 | sed -e "s/^"$DIRPREFIX"0*//"`
AUX=$((AUX+1))
AUX=`printf "%06d\n" $AUX`
SIMDIR="simulations/"$DIRPREFIX$AUX
NTHREADS=
COMPTYPE=FAST_OMP	# See Makefile:
			# STD: Standar compilation
			# FAST: With Ofast optimizations
			# OPT: With O3 optimizations
			# DBG: With debug options
			# OMP: With OpenMP
			# FAST_OMP: With Ofast optimizations and OpenMP
			# O3_OMP: With O3 optimizations and OpenMP

while getopts "d:n:c:" o; do
	case "${o}" in
		d)
			SIMDIR=$OPTARG
			;;
		n)
			NTHREADS=$OPTARG
			;;
		c)
			COMPTYPE=$OPTARG
			;;
		\?) 
			echo "Invalid option -$OPTARG" >&2
			break
			;;
	esac
done

if [ -d $SIMDIR ]; then
	echo "Error: directory already exists"
	exit 1
elif [ $SIMDIR == "" ]; then
	echo "Error: you must define a directory"
	exit 1
fi

mkdir -p $SIMDIR
mkdir -p $SIMDIR/data
mkdir -p $SIMDIR/src
mkdir -p $SIMDIR/bin
cp src/*.cpp src/*.h src/Makefile $SIMDIR/src/
make COMP_TYPE=$COMPTYPE -C $SIMDIR/src
mv $SIMDIR/src/$BIN $SIMDIR/bin/$BIN-$AUX

if [ "$NTHREADS" != "" ]; then
	export OMP_NUM_THREADS=$NTHREADS
fi
ulimit -s unlimited
cd $SIMDIR/data
time ../bin/$BIN-$AUX > out.dat &


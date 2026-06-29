#!/bin/bash

# CONTROL
EXECUTABLE="./build/NSSolver"
INPUT_DIR="./ioSrc"

# CORES
CORES=${1:-6}
if [ "$CORES" -gt 10 ]; then
	echo "Core limit exceeded. Capping at 8."
	CORES=8
fi

# RANGE
if [ -z "$2" ]; then
	SIMRANGE="*"
else
	ARGS_RANGE=$(echo "$@" | cut -d' ' -f2- | tr ' ' ',')
	SIMRANGE="{$ARGS_RANGE}"
fi

# HEADER
rm -f simulation_log.txt
echo "----------------------------------------------------" > simulation_log.txt
echo "Running simulations w/ $CORES cores ... $(date)" | tee -a simulation_log.txt
echo "----------------------------------------------------" >> simulation_log.txt

# RUN SIMULATIONS
eval ls "$INPUT_DIR"/exNS"$SIMRANGE".json | sort -V | parallel -u -j "$CORES" \
	"$EXECUTABLE {} && echo Finished {} at \$(date '+%H:%M:%S') | tee -a simulation_log.txt"

# FOOTNOTE
echo "----------------------------------------------------" >> simulation_log.txt
echo "All simulations complete ... $(date)" | tee -a simulation_log.txt



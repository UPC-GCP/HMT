#!/bin/bash

### CONTROL
EXECUTABLE="./build/ROCK"
INPUT_DIR="./ioSrc"
INPUT_FILE=${1}

### CORES
CORES=${2:-6}
if [ "$CORES" -gt 10 ]; then
	echo "Core limit exceeded. Capping at 8."
	CORES=8
fi

### HEADER
rm -f simulation_log.txt
echo "----------------------------------------------------" > simulation_log.txt
echo "Running simulations w/ $CORES cores ... $(date)" | tee -a simulation_log.txt
echo "----------------------------------------------------" >> simulation_log.txt

### RUN SIMULATIONS
eval ls "$INPUT_DIR"/"$INPUT_FILE"*.json | sort -V | parallel -u -j "$CORES" \
	"$EXECUTABLE {} && echo Finished {} at \$(date '+%H:%M:%S') | tee -a simulation_log.txt"

### FOOTNOTE
echo "----------------------------------------------------" >> simulation_log.txt
echo "All simulations complete ... $(date)" | tee -a simulation_log.txt

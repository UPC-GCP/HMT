#!/bin/bash

### CONTROL
### GENERAL VARIABLE DEFINITION
EXECUTABLE="./build/ROCK"
INPUT_DIR="./ioSrc"
BASE_FILE="exSpectral.json" # ex[SOLVERTYPE]
OUT_BASE="cfgBurgers" # cfg[CASE]
# BASE FILE MODIFIED INTO OUT_BASE[PARAMETERS]
# sample files: ex[SOLVERTYPE], configuration files: cfg[CASE][PARAMETERS]

### RANGES
### ALL CHANGHING VALUES LISTED HERE (SWEEPVARIABLE=("Opt 1", "Opt 2", ... ))
SWEEPSOLVER=("DNS" "LES")
SWEEPN=(20 40 60 80 100) # Only for DNS
SWEEPCS=(0.4223 0.35 0.25 0.15 0.05) # Only for LES

# Parameter loop
for solver in "${SWEEPSOLVER[@]}"; do # SOLVER loop
	
    # DNS Only (SWEEP N)
    if [[ $solver == "DNS" ]]; then
        for N in "${SWEEPN[@]}"; do
            # Name
            FILENAME="$INPUT_DIR/${OUT_BASE}${solver}N${N}.json"
            probes="[
                {\"type\": \"Modal\", \"t\": [0, 100], \"x0\": [0], \"x1\": [$N], \"nWrite\": null}
            ]"
            echo "${OUT_BASE}${solver}N${N}.json"


            # Save
            jq --arg argSol "$solver" \
                --argjson argN "$N" \
                --argjson argProbe "$probes" \
                '.solver = $argSol | .N = $argN | .probes = $argProbe' \
                "$INPUT_DIR/$BASE_FILE" > "$FILENAME"
        done
    fi

    # LES Only (SWEEP CS)
    if [[ $solver == "LES" ]]; then
        for Cs in "${SWEEPCS[@]}"; do
            # File
            FILENAME="$INPUT_DIR/${OUT_BASE}${solver}Cs${Cs}.json"
            echo "${OUT_BASE}${solver}Cs${Cs}.json"

            # Save
            jq --arg argSol "$solver" \
                --argjson argCs "$Cs" \
                '.solver = $argSol | .Cs = $argCs' \
                "$INPUT_DIR/$BASE_FILE" > "$FILENAME"
        done
    fi

done

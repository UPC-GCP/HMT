#!/bin/bash

### CONTROL
### GENERAL VARIABLE DEFINITION
EXECUTABLE="./build/ROCK"
INPUT_DIR="./ioSrc"
BASE_FILE="exSpectral.json" # ex[SOLVERTYPE]
OUT_BASE="cfgBurgers" # cfg[CASE]
# BASE FILE MODIFIED INTO OUT_BASE[PARAMETERS]
# sample files: ex[SOLVERTYPE], configuration files: cfg[CASE][PARAMETERS]
# [PARAMETERS]: Include the variables analyzed in the sweep [SOLVER+VAR+VALUE]

### RANGES
### ALL CHANGHING VALUES LISTED HERE (SWEEPVARIABLE=("Opt 1", "Opt 2", ... ))
SWEEPSOLVER=("DNS" "LES")
SWEEPDT=(0.1 0.01 0.001 0.0001) # DNS: N = 100, LES: N = 20
SWEEPC1=(0.1 0.45 0.8 1.15 1.5) # DNS: N = 100, LES: N = 20

# DNS
SWEEPDNSN=(20 40 60 80 100) # dt = 0.0001

#LES
SWEEPLESN=(10 15 20 25 30) # dt = 0.0001
SWEEPCK=(0.05 0.25 0.40 0.4223 0.45 0.60 0.80) # N = 20, dt = 0.0001

# Parameter loop
for solver in "${SWEEPSOLVER[@]}"; do # SOLVER loop
	
    # DNS
    if [[ $solver == "DNS" ]]; then
        # Temporal Convergence - T
        for dt in "${SWEEPDT[@]}"; do
            # Name
            FILENAME="$INPUT_DIR/${OUT_BASE}TempT${solver}dt${dt}.json"
            echo "$FILENAME"

            # Parameters
            N=100
            probes="[
                {\"type\": \"Modal\", \"t\": [0, 100], \"x0\": [0], \"x1\": [$N], \"nWrite\": null}
            ]"

            # Create .json
            jq --arg argSol "$solver" \
                --argjson argDT "$dt" \
                --argjson argN "$N" \
                --argjson argProbe "$probes" \
                '.solver = $argSol | .timeStep = $argDT | .N = $argN | .probes = $argProbe' \
                "$INPUT_DIR/$BASE_FILE" > "$FILENAME"
        done

        # Temporal Convergence - C1
        for C1 in "${SWEEPC1[@]}"; do
            # Name
            FILENAME="$INPUT_DIR/${OUT_BASE}TempC${solver}C${C1}.json"
            echo "$FILENAME"

            # Parameters
            dt=1
            N=100
            probes="[
                {\"type\": \"Modal\", \"t\": [0, 100], \"x0\": [0], \"x1\": [$N], \"nWrite\": null}
            ]"

            # Create .json
            jq --arg argSol "$solver" \
                --argjson argC "$C1" \
                --argjson argDT "$dt" \
                --argjson argN "$N" \
                --argjson argProbe "$probes" \
                '.solver = $argSol | .C1 = $argC | .timeStep = $argDT | .N = $argN | .probes = $argProbe' \
                "$INPUT_DIR/$BASE_FILE" > "$FILENAME"
        done

        # Mesh Refinement - N [READY]
        for N in "${SWEEPDNSN[@]}"; do
            # Name
            FILENAME="$INPUT_DIR/${OUT_BASE}Mesh${solver}N${N}.json"
            echo "$FILENAME"

            # Parameters
            probes="[
                {\"type\": \"Modal\", \"t\": [0, 100], \"x0\": [0], \"x1\": [$N], \"nWrite\": null}
            ]"
            dt=0.0001

            # Create .json 
            jq --arg argSol "$solver" \
                --argjson argN "$N" \
                --argjson argProbe "$probes" \
                --argjson argDT "$dt" \
                '.solver = $argSol | .N = $argN | .probes = $argProbe | .timeStep = $argDT' \
                "$INPUT_DIR/$BASE_FILE" > "$FILENAME"
        done
    fi

    # LES
    if [[ $solver == "LES" ]]; then
        # Temporal Convergence - T [READY]
        for dt in "${SWEEPDT[@]}"; do
            # Name
            FILENAME="$INPUT_DIR/${OUT_BASE}TempT${solver}dt${dt}.json"
            echo "$FILENAME"

            # Parameters
            N=20
            probes="[
                {\"type\": \"Modal\", \"t\": [0, 100], \"x0\": [0], \"x1\": [$N], \"nWrite\": null}
            ]"

            # Create .json
            jq --arg argSol "$solver" \
                --argjson argDT "$dt" \
                --argjson argN "$N" \
                --argjson argProbe "$probes" \
                '.solver = $argSol | .timeStep = $argDT | .N = $argN | .probes = $argProbe' \
                "$INPUT_DIR/$BASE_FILE" > "$FILENAME"
        done

        # Temporal Convergence - C1 [READY]
        for C1 in "${SWEEPC1[@]}"; do
            # Name
            FILENAME="$INPUT_DIR/${OUT_BASE}TempC${solver}C${C1}.json"
            echo "$FILENAME"

            # Parameters
            dt=1
            N=20
            probes="[
                {\"type\": \"Modal\", \"t\": [0, 100], \"x0\": [0], \"x1\": [$N], \"nWrite\": null}
            ]"

            # Create .json
            jq --arg argSol "$solver" \
                --argjson argC "$C1" \
                --argjson argDT "$dt" \
                --argjson argN "$N" \
                --argjson argProbe "$probes" \
                '.solver = $argSol | .C1 = $argC | .timeStep = $argDT | .N = $argN | .probes = $argProbe' \
                "$INPUT_DIR/$BASE_FILE" > "$FILENAME"
        done

        # Mesh Refinement - N [READY]
        for N in "${SWEEPLESN[@]}"; do
            # Name
            FILENAME="$INPUT_DIR/${OUT_BASE}Mesh${solver}N${N}.json"
            echo "$FILENAME"

            # Parameters
            probes="[
                {\"type\": \"Modal\", \"t\": [0, 100], \"x0\": [0], \"x1\": [$N], \"nWrite\": null}
            ]"
            dt=0.0001

            # Create .json 
            jq --arg argSol "$solver" \
                --argjson argN "$N" \
                --argjson argProbe "$probes" \
                --argjson argDT "$dt" \
                '.solver = $argSol | .N = $argN | .probes = $argProbe | .timeStep = $argDT' \
                "$INPUT_DIR/$BASE_FILE" > "$FILENAME"
        done

        # Kolmogorov Constant [READY]
        for CK in "${SWEEPCK[@]}"; do
            # File
            FILENAME="$INPUT_DIR/${OUT_BASE}KolC${solver}CK${CK}.json"
            echo "$FILENAME"

            # Parameters
            dt=0.0001
            N=20
            probes="[
                {\"type\": \"Modal\", \"t\": [0, 100], \"x0\": [0], \"x1\": [$N], \"nWrite\": null}
            ]"

            # Create .json 
            jq --arg argSol "$solver" \
                --argjson argCK "$CK" \
                --argjson argDT "$dt" \
                --argjson argN "$N" \
                --argjson argProbe "$probes" \
                '.solver = $argSol | .CK = $argCK | .timeStep = $argDT | .N = $argN | .probes = $argProbe' \
                "$INPUT_DIR/$BASE_FILE" > "$FILENAME"
        done
    fi

done

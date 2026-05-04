#!/bin/bash

### CONTROL
EXECUTABLE="./build/MainSolver"
INPUT_DIR="./RunConfig"
BASE_FILE="4Materials.json"
counter=0


### RANGES
sSWEEP=("implicit" "crank-nicolson" "explicit")
tSWEEP=(0.5 1 2 5)
NSWEEP=(50 100 200 400)


# scheme loop
for scheme in "${sSWEEP[@]}"; do
	
	# timeStep Loop
	for dt in "${tSWEEP[@]}"; do
	
		# Control
		((counter++))
		echo "Case ${counter}: ${scheme}_t${dt}"

		# Save .json
		FILENAME="$INPUT_DIR/Case_${counter}.json"
		jq --arg s "$scheme" \
			--argjson t "$dt" \
			'.scheme = $s | .timeStep = $t' \
			"$INPUT_DIR/$BASE_FILE" > "$FILENAME"
	
	done

	# N Loop
	for N in "${NSWEEP[@]}"; do

		# Control
		((counter++))
		echo "Case ${counter}: ${scheme}_N${N}"

		# Values
		rN0=$(((5 * N + 5) / 11))
		# rN1=$(((6 * N + 5) / 11))
		rN1=$((N - rN0))
		rN2=$(((4 * N + 4) / 8))
		rN3=$(((3 * N + 4) / 8))
		# rN4=$(((1 * N + 4) / 8))
		rN4=$((N - rN2 - rN3))

		refinement="[
			{\"axis\": 0, \"N\": $rN0, \"range\": [0, 0.5]},
			{\"axis\": 0, \"N\": $rN1, \"range\": [0.5, 1.1]},
			{\"axis\": 1, \"N\": $rN2, \"range\": [0, 0.4]},
			{\"axis\": 1, \"N\": $rN3, \"range\": [0.4, 0.7]},
			{\"axis\": 1, \"N\": $rN4, \"range\": [0.7, 0.8]}
		]"
		
		gridShape="[$N, $N]"

		# Save .json
		FILENAME="$INPUT_DIR/Case_${counter}.json"
		jq --arg s "$scheme" \
			--argjson shape "$gridShape" \
			--argjson ref "$refinement" \
			'.scheme = $s | .N = $shape | .refinement = $ref' \
			"$INPUT_DIR/$BASE_FILE" > "$FILENAME"

	done
	
done


# CLEAN NAMES
for f in "$INPUT_DIR"/Case_[0-9].json; do mv "$f" "${f/Case_/Case_0}"; done


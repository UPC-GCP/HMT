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

# BUILD CHECK
#echo "--- Compiling Project ---"
#cmake --build build --config Release
#if [ $? -ne 0 ]; then
#    echo "Build failed! Exiting."
#    exit 1
#fi

# CREATE FILES
for scheme in "${sSWEEP[@]}"; do
	
	# timeStep Loop
	for dt in "${tSWEEP[@]}"; do
	
		# Control
		echo "Case: ${scheme}_t${dt}"
		((counter++))

		# Save .json
		FILENAME="$INPUT_DIR/Case_${counter}.json"
		jq --arg s "$scheme" \
			--argjson t "$dt" \
			'.scheme = $s | .timeStep = $t' \
			"$INPUT_DIR/$BASE_FILE" > "$FILENAMEn"
	
	done

	# N Loop
	for N in "${NSWEEP[@]}"; do

		# Control
		echo "Case: ${scheme}_N${N}"
		((counter++))

		# Values
		rN0=$(((5 * N + 5) / 11))
		rN1=$(((6 * N + 5) / 11))
		rN2=$(((4 * N + 4) / 8))
		rN3=$(((3 * N + 4) / 8))
		rN4=$(((1 * N + 4) / 8))

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
			"$BASE_FILE" > "$FILENAME"


	done
	
done


#!/bin/bash

# CONTROL
iFiles=${1:-0}

# OPEN FILES
if [ "$iFiles" -eq 0 ]; then
	vim -p ./RunConfig/example.json MainSolver.cpp Mesh.cpp Mesh.h Discretizer.cpp Discretizer.head
elif [ "$iFiles" -eq 1 ]; then
	vim -p MainSolver.cpp ./RunConfig/example.json Mesh.cpp Mesh.h ExpressionParser.cpp ExpressionParser.h
else
	echo "No file group recognized. (Max: 1)"
fi


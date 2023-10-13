#!/bin/bash

gcc timing.c -I./include ./src/matrix.c -o time_matvec
if [ ! -d "./docs" ]; then
	mkdir ./docs
fi
./time_matvec 1000 1000 10 > docs/timing.txt
./time_matvec 2000 2000 10 >> docs/timing.txt
./time_matvec 3000 3000 10 >> docs/timing.txt
./time_matvec 4000 4000 10 >> docs/timing.txt

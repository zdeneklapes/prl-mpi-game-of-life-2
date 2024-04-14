#!/bin/bash

game_input=$1
steps=$2

function build_mpi() {
    program="life"
    mpic++ -g -O0 --prefix /usr/local/share/OpenMPI -o "$program" "$program".cpp || die "mpic++ failed"
}

function run_mpi() {
    proc=3
    build_mpi || die "build failed"
    program="life"
#    echo "Pocet procesoru: $proc; pocet kroku: $steps"
#    echo "================"
    mpirun -np "$proc" ./"$program" "$game_input" "$steps" || die "mpirun failed"
}

build_mpi
run_mpi

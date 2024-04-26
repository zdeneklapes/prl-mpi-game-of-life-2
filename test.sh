#!/bin/bash

# Wrap-around version of Game of Life

game_input=$1
steps=$2

RED='\033[0;31m'
NC='\033[0m'

function die() {
    # Print error message on stdout and exit
    printf "${RED}ERROR: $1${NC}\n"
    exit 1
}

function build_mpi() {
    program="life"
    mpic++ -g -O0 --prefix /usr/local/share/OpenMPI -o "$program" "$program".cpp || die "mpic++ failed"
}

function run_mpi() {
    rows_in_file=$(wc -l < "$game_input")
    available_proceses=$(nproc)
    if [ "$rows_in_file" -gt "$available_proceses" ]; then
        proc="$available_proceses"
    else
        proc="$rows_in_file"
    fi
#    proc=6 # accept any number of processes available on the system
    build_mpi || die "build failed"
    program="life"

#    echo "Pocet procesoru: $proc; pocet kroku: $steps"
#    echo "================"
    mpirun -np "$proc" ./"$program" "$game_input" "$steps" || die "mpirun failed"
}

build_mpi
run_mpi

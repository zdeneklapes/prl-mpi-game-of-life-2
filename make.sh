#!/bin/bash

###########################################################
###########################################################
RM="rm -rfd"
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'

###########################################################
###########################################################
function push_environs() {
    # Push environment variables to the server
    CMD="npx dotenv-vault@latest push"
    for stage in "development" "production" "ci" "staging"; do
        CMD_FINAL="${CMD} ${stage}"
        if [ $DEBUG -eq 1 ]; then echo "$CMD_FINAL"; else eval "$CMD_FINAL"; fi
    done
}

function delete_local_environs() {
    # Delete environment variables from the server
    CMD="rm"
    for stage in ".development" ".production" ".ci" ".staging" ""; do
        CMD_FINAL="${CMD} .env${stage}"
        if [ $DEBUG -eq 1 ]; then echo "$CMD_FINAL"; else eval "$CMD_FINAL"; fi
    done
}

function pull_environs() {
    # Pull environment variables from the server
    CMD="npx dotenv-vault@latest pull"
    for stage in "development" "production" "ci" "staging"; do
        CMD_FINAL="${CMD} ${stage}"
        if [ $DEBUG -eq 1 ]; then echo "$CMD_FINAL"; else eval "$CMD_FINAL"; fi
    done
}

###########################################################
###########################################################
function clean() {
    # Clean project folder in order to see what will be done, set env variable $DEBUG=1
    ${RM} *.zip
    # Folders
    for folder in \
            "venv" \
            "*__pycache__" \
            "*.ruff_cache" \
            "*.pytest_cache" \
            "*.cache" \
            "*htmlcov*" \
            "skip-covered"\
            ; do
        if [ "$DEBUG" -eq 1 ]; then find . -type d -iname "${folder}"; else find . -type d -iname "${folder}" | xargs ${RM} -rf; fi
    done
    # Files
    for file in \
            "*.DS_Store" \
            "tags" \
            "db.sqlite3" \
            "*.png" \
            "*.zip" \
            "*.log" \
            "coverage.xml" \
            "*.coverage" \
            "coverage.lcov" \
            ; do
        if [ "$DEBUG" -eq 1 ]; then find . -type f -iname "${file}"; else find . -type f -iname "${file}" | xargs ${RM}; fi
    done
}

function build() {
    mkdir -p build || die "mkdir failed"
    cd build || die "cd failed"
    cmake .. || die "cmake failed"
    make || die "make failed"
    cd .. || die "cd failed"
}

function build_mpi() {
    program="life"
    mpic++ -g -O0 --prefix /usr/local/share/OpenMPI -o "$program" "$program".cpp || die "mpic++ failed"
}

function run_mpi() {
    steps=1
    proc=3
    build_mpi || die "build failed"
    program="life"
    echo "Pocet procesoru: $proc; pocet kroku: $steps"
    echo "================"
    mpirun -np "$proc" ./"$program" ./tests/2-input-1-steps.txt "$steps" || die "mpirun failed"
}


function valgrind_mpi() {
    build_mpi || die "build_mpi failed"
    valgrind --leak-check=yes mpirun -n 2 ./life ./tests/1-input.txt 1
}

function rsync_to_server() {
    # Rsync to the server
    # ENVIRONMENT VARIABLES
    #   - PROJECT_PATH
    #   - SERVER_URL
    #   - USER
    #   - DEBUG

    # Set environment variables
    if [ -z "$PROJECT_PATH" ]; then PROJECT_PATH="~/repos/prl-1"; fi
#    if [ -z "$SERVER_URL" ]; then SERVER_URL="merlin.fit.vutbr.cz"; fi # eva.fit.vutbr.cz
    if [ -z "$SERVER_URL" ]; then SERVER_URL="eva.fit.vutbr.cz"; fi
    if [ -z "$USER_SERVER" ]; then USER_SERVER="xlapes02"; fi

    CMD1="mkdir -p '${PROJECT_PATH}'"

    ssh "${USER_SERVER}@${SERVER_URL}" "${CMD1}"

    rsync --archive --verbose --compress --delete \
        pms.cpp \
        CMakelists.txt \
        docker-compose.yml Dockerfile \
        make.sh README.md generateNums.sh \
        "${USER_SERVER}@${SERVER_URL}:${PROJECT_PATH}"
}

function help() {
    # Print usage on stdout
    echo "Available functions:"
    for file in "make.sh"; do
        function_names=$(cat ${file} | grep -E "(\ *)function\ +.*\(\)\ *\{" | sed -E "s/\ *function\ +//" | sed -E "s/\ *\(\)\ *\{\ *//")
        for func_name in ${function_names[@]}; do
            printf "    $func_name\n"
            awk "/function ${func_name}()/ { flag = 1 }; flag && /^\ +#/ { print \"        \" \$0 }; flag && !/^\ +#/ && !/function ${func_name}()/  { print "\n"; exit }" ${file}
        done
    done
}

function usage() {
    # Print usage on stdout
    help
}

function die() {
    # Print error message on stdout and exit
    printf "${RED}ERROR: $1${NC}\n"
    exit 1
}

function main() {
    # Main function: Call other functions based on input arguments
    [[ "$#" -eq 0 ]] && die "No arguments provided"
    while [ "$#" -gt 0 ]; do
        "$1" || die "Unknown function: $1()"
        shift
    done
}
main "$@"

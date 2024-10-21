#!/bin/bash

# Use this script to install the example (i.e. create the executable)
# Use argument 0 to re-build programs and libraries that have been modified sice last build
# Use argument 1 for complete re-build (removing everythinh from build dir and start over)

ERROR_FLAG=0

# Handeling the argument passed to the script
MODE=$1
if [ -z "$1" ]; then
    echo "Useage " $0 " 0 <re-build modified bin / libs> 1 <re-build all>"
    ERROR_FLAG=1
fi
if [ $ERROR_FLAG == 1 ]; then
    exit 1
fi

# This forces a complete re-build
if [ $MODE == 1 ]; then
    rm -r ./build/*
fi

# Now we run cmake
cd build
#cmake -DGeant4_DIR=/home/sdekockere/geant4/geant4.10.05-install/lib/Geant4-10.5.1 ../
cmake ../

# Now we make
make

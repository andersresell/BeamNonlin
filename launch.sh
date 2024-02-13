#!/bin/bash

clear
source build.sh

if [ "$BUILD_TYPE" == "release" ]; then
    builddir=$BeamNonlinHome/build_release
elif [ "$BUILD_TYPE" == "debug" ]; then
    builddir=$BeamNonlinHome/build_debug
else
    echo "error: invalid build type"
    exit 1
fi

if [ $? -ne 0 ]; then
    exit 1
fi


simdir=$BeamNonlinHome/testing
cd $simdir
$builddir/BeamNonlin 


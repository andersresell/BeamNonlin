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
$builddir/BeamNonlin cantilever.yml #> $BeamNonlinHome/console_output.txt

# echo "Writing profiling output.."
# gprof $builddir/BeamNonlin gmon.out | gprof2dot -w -s | dot -Gdpi=200 -Tpng -o profiling_output.png
# rm gmon.out
# echo "Done"

#valgrind --tool=memcheck $builddir/BeamNonlin cantilever.yml
#valgrind --tool=callgrind --callgrind-out-file=callgrind.out $builddir/BeamNonlin cantilever.yml
#kcachegrind callgrind.out
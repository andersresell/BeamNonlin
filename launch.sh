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

input_file=cantilever-2d-steady.yml
input_file=cantilever-2d-small.yml
input_file=cantilever-small-problem-unstable.yml
input_file=spinning-top.yml
input_file=cantilever.yml



$builddir/BeamNonlin $input_file #> $BeamNonlinHome/console_output.txt

# echo "Writing profiling output.."
# gprof $builddir/BeamNonlin gmon.out | gprof2dot -w -s | dot -Gdpi=200 -Tpng -o profiling_output.png
# rm gmon.out
# echo "Done"

#valgrind --tool=memcheck $builddir/BeamNonlin cantilever.yml
#valgrind --tool=callgrind --callgrind-out-file=callgrind.out $builddir/BeamNonlin cantilever.yml
#kcachegrind callgrind.out
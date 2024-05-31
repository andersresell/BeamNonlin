#!/bin/bash
clear

run_plotter=0
input_file_specififed=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        -p|--plot)
            run_plotter=1
            echo "Plotter will run after simulation is complete"
            ;;
        -f=*|--file=*)
            input_file="${1#*=}"
            input_file_specififed=1
            ;;
        *)
            echo "Invalid option"
            exit 1
            ;;
    esac
    shift
done

if [ $input_file_specififed -eq 0 ]; then
    echo "Input file has not been specified, use -f=<file> or --file=<file>"
    exit 1
fi

#Sourcing the build script with no command line options, which should result in release build
source $BeamNonlinHome/build.sh 

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


#simdir=$BeamNonlinHome/testing 
#cd $simdir
simdir=$(pwd)
echo "Simulation output is saved in "$simdir

$builddir/BeamNonlin $input_file #> $BeamNonlinHome/console_output.txt

if [ $? -eq 0 ]; then
    echo "simulation returned without errors"
else
    echo "simulation returned with errors"
    exit 1
fi

if [ $run_plotter -eq 1 ]; then
    python $BeamNonlinHome/plotter.py
fi






# echo "Writing profiling output.."
# gprof $builddir/BeamNonlin gmon.out | gprof2dot -w -s | dot -Gdpi=200 -Tpng -o profiling_output.png
# rm gmon.out
# echo "Done"

#valgrind --tool=memcheck $builddir/BeamNonlin cantilever.yml
#valgrind --tool=callgrind --callgrind-out-file=callgrind.out $builddir/BeamNonlin cantilever.yml
#kcachegrind callgrind.out
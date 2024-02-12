#!/bin/bash
clear

rm -f $DrillSimHome/build_debug/BeamNonlin
rm -f $DrillSimHome/build_release/BeamNonlin


while [ "$#" -gt 0 ]; do
    case "$1" in
        -g)
            export BUILD_TYPE="debug"
            ;;
        -c)
            make -C $DrillSimHome/src clean
            exit $?
            ;;
        *)
            echo "Invalid build option: $1"
            exit 1
            ;;
    esac
    shift
done

make -C $DrillSimHome/src all

if [ $? -ne 0 ]; then
    exit 1
fi

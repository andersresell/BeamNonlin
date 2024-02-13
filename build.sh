#!/bin/bash

#not sure why it isnt sufficient to specify this in .bashrc, when using it in vscode tasks
export BeamNonlinHome="/home/anders/projects/BeamNonlin" 

rm -f $BeamNonlinHome/build_debug/BeamNonlin
rm -f $BeamNonlinHome/build_release/BeamNonlin

export BUILD_TYPE="release"


while [ "$#" -gt 0 ]; do
    case "$1" in
        -g)
            export BUILD_TYPE="debug"
            ;;
        -c)
            make -C $BeamNonlinHome/src clean
            exit $?
            ;;
        *)
            echo "Invalid build option: $1"
            exit 1
            ;;
    esac
    shift
done

echo "HOME $BeamNonlinHome"
make -C $BeamNonlinHome/src all

if [ $? -ne 0 ]; then
    exit 1
fi

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
        -c|--clean)
            make -C $BeamNonlinHome/src clean
            exit $?
            ;;
        uf|--user-force=*)
            user_force_function_file="${1#*=}"
            export USER_DEFINED_EXTERNAL_FORCE_FILE="$(pwd)/$user_force_function_file"
            ;;
        *)
            echo "Invalid build option: $1"
            exit 1
            ;;
    esac
    shift
done

make -C $BeamNonlinHome/src all

if [ $? -ne 0 ]; then
    exit 1
else
    echo "Build successful"
fi

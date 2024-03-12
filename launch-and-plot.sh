#!/bin/bash
clear

./launch.sh 

if [ $? -eq 0 ]; then
    echo "simulation returned without errors"
else
    echo "simulation returned with errors"
    exit 1
fi

 
python plotter.py
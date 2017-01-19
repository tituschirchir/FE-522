#!/usr/bin/env bash
if [ -d "build" ]; then
    echo "Build directory already exists! Do you want to recompile? [Y/N]";   read response
    if [ "$response" = "y" -o "$response" = "yes" -o  "$response" = "Y" -o "$response" = "Yes" ]; then
        rm -rf build;
        mkdir build;
        cd build;
        cmake ../;
        make;
        cd ..;
    fi
else
    mkdir build;
    cd build;
    cmake ../;
    make;
    cd ..;
fi
./build/brownianmotion;

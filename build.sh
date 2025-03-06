#!/bin/bash

cmake -S . -B build
cmake --build build/
cp build/libGGFW* .
root -l exec/test.C
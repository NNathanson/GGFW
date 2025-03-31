#!/bin/bash

cmake -S . -B build
cmake --build build/
cp build/libGGFW* .
root -l -q exec/test.C
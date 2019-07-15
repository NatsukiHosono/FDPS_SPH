#!/bin/sh
cd build
make 
cd ..
./sph.out -i input/input.txt

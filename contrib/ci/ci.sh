#!/bin/bash

bash ./CloneFDPS.sh
mkdir build
cd build
cmake ..
make

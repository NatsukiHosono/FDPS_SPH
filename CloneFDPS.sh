#!/bin/sh
git clone https://github.com/FDPS/FDPS.git
cd FDPS

# This checks out a version of FDPS that is known to work with FDPS_SPH
# Modify if you want to use a different FDPS version
git checkout v5.0c
cd ..

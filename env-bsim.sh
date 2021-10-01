#!/bin/bash

# Shell script for setting up BSim environment on Graham compute cluster.
# Newest available version of java (13.0.2) is loaded in place of default (1.8.0).
# Creates alias for building bsim by executing 'bsim-build' in command line.

echo "Setting up BSim environment..."
module load ant
module load openmpi
module load mpi4py #needed for multiprocessing in python and abc code
module load StdEnv/2020 
module load scipy-stack/2020b
module load opencv

source $HPC_BSIM/python-env.sh


export BSIM_DIR=$PWD
alias bsim-build='ant -f $BSIM_HOME/bsim-build-tree.xml'

echo "Done."

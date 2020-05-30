#!/bin/sh

###############################################################################
# Copyright (c) The JETSCAPE Collaboration, 2018
#
# For the list of contributors see AUTHORS.
#
# Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
#
# or via email to bugs.jetscape@gmail.com
#
# Distributed under the GNU General Public License 3.0 (GPLv3 or later).
# See COPYING for details.
##############################################################################

if [ -z $LD_LIBRARY_PATH ]; then
export LD_LIBRARY_PATH
fi

if [ -z $DYLD_LIBRARY_PATH ]; then 
export DYLD_LIBRARY_PATH
fi

export BASEDIR=${HOME}
export PYTHIAINSTALLDIR=/home/td115/lib/pythia8235

export JetScape=${PWD}/lib
export LD_LIBRARY_PATH=${JetScape}:${LD_LIBRARY_PATH}
#only for Mac needed
export DYLD_LIBRARY_PATH=${JetScape}:${DYLD_LIBRARY_PATH}

# PYTHIA8 directory
export PYTHIA8DIR=${PYTHIAINSTALLDIR}
export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}
#export PYTHIA8_INCLUDE_DIR=`${PYTHIA8DIR}/bin/pythia8-config --includedir`/Pythia8
#export PYTHIA8_LIBRARIES=`${PYTHIA8DIR}/bin/pythia8-config --libdir`
export LD_LIBRARY_PATH=${PYTHIA8DIR}/lib:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${PYTHIA8DIR}/lib:${DYLD_LIBRARY_PATH}

#ROOT setup
export ROOTSYS="/opt/local/libexec/root5"
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}

if [ -z ${TERM} -o -z ${SHELL} ]; then exit 0
fi

export EIGEN_INSTALL_DIR=/home/td115/lib/eigen-eigen-3.3.7
export EIGEN3_ROOT=/home/td115/lib/eigen-eigen-3.3.7
export GSL=$(gsl-config --prefix)
export GSL_HOME=$(gsl-config --prefix)
export GSL_ROOT_DIR=$(gsl-config --prefix)
export JETSCAPE_DIR=`readlink -f .`
export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code
export number_of_cores=`nproc --all`
export CC=gcc
export CXX=g++
export OpenMP_CXX=g++

echo ''
echo 'Setup JetScape Library'
echo '======================'
echo ''
echo "<I>---------------Info--------------------<I>"
echo "Setting up the following environments: "
echo "JetScape: " ${JetScape}
echo "Pythia8: "${PYTHIA8DIR}/lib    
echo "<I>---------------Info--------------------<I>"
echo ""

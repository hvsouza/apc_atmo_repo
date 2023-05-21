#!/bin/bash

DUNESW_version=v09_72_00d00

QUALS=e20:prof
DIRECTORY=larsoft
USERNAME=`whoami`
export WORKDIR=/dune/app/users/${USERNAME}/atmo_analysis
if [ ! -d "$WORKDIR" ]; then
  export WORKDIR=`echo ~`
  echo $WORKDIR
  exit 0
fi

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunesw ${DUNESW_version} -q ${QUALS}

cd ${WORKDIR}
touch ${DIRECTORY}
rm -rf ${DIRECTORY}
mkdir ${DIRECTORY}
cd ${DIRECTORY}
mrb newDev -q ${QUALS}
source ${WORKDIR}/${DIRECTORY}/localProducts*/setup
mkdir work
cd srcs

mrb g -t ${DUNESW_version} duneana

cd ${MRB_BUILDDIR}
mrbsetenv
mrb i -j16 --generator ninja
# then just `ninja install` for faster build+compile

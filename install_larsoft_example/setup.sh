DIRECTORY=larsoft
USERNAME=`whoami`

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
export WORKDIR=/dune/app/users/${USERNAME}/analysis_pierre

if [ ! -d "$WORKDIR" ]; then
  export WORKDIR=`echo ~`
fi

DUNESW_version=v09_72_00d00
QUALS=e20:prof
setup dunesw ${DUNESW_version} -q ${QUALS}

cd $WORKDIR/$DIRECTORY
source localProducts*/setup
cd work
mrbslp
mrbsetenv

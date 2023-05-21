# Example on installing larsoft


## Simple and easy

If you are just joing to run LArSoft (without modifying anything). You can just run the following commands:

``` sh
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
DUNESW_version=v09_72_00d00
QUALS=e20:prof
setup dunesw ${DUNESW_version} -q ${QUALS}
```

And you are good to be happy!

## Simple but painful 

In case you want to modify modules or services in larsoft, this is how you do it:

This example will install larsoft in the folder `atmo_analysis/folder`. The module `duneana` is installed built together, this is done in case you want to modify something in larsoft. Otherwise you don't need this.

Create the folder:

``` sh
USERNAME=`whoami`
mkdir /dune/app/users/${USERNAME}/atmo_analysis
```

Before running the script `newDev.sh`, **PLEASE** ssh in the building machines: `ssh dunebuild01.fnal.gov`.

Now you can run `source newDev.sh`. After installing, `exit` from the building machines

If everything works fine, you need to `source` the script `setup.sh` every time you have a new terminal open.

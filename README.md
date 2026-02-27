# UNICICLES

## Build notes

### Ubuntu 22.04 (serial, debug)

Assumes BISICLES at $BISICLES_HOME/bisicles-uob and Chombo at $BISICLES_HOME/Chombo.
Both built with DEBUG=TRUE, OPT=FALSE, MPI=FALSE, USE_PETSC=FALSE

```

cd $BISICLES_HOME
git clone git@github.com:ggslc/unicicles.git
cd unicicles

# Build glimmer-cism
cd glimmer-cism
./bootstrap
cd ..
mkdir glimbike-serial
cd glimbike-serial

BIKE=$BISICLES_HOME/bisicles-uob/
BIKE_CONFIG=2d.Linux.64.g++.gfortran.DEBUG
HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/
NETCDF=/usr

PYTHON=python3 FC=gfortran FCFLAGS="-fno-range-check -ffree-line-length-0 -DBISICLES_CDRIVER -DNO_RESCALE -g -I$BIKE/code/src " LDFLAGS="-L$BIKE/code/lib -lBisicles$BIKE_CONFIG -lChomboLibs$BIKE_CONFIG -lpython3.10 -L$HDF5 -lhdf5 -lz " ../glimmer-cism/configure --with-netcdf=$NETCDF --with-hdf5=$HDF5 --prefix=$PWD --disable-python
make
make install

# build wrappers/ukesm-ice_NETCDF
cd wrappers/ukesm-ice_NETCDF
make clean -f Makefile.ubuntu22.04
make -f Makefile.ubuntu22.04

```

### Ubuntu 22.04 (parallel, opt, petsc, debug)

Assumes BISICLES at $BISICLES_HOME/bisicles-uob and Chombo at $BISICLES_HOME/Chombo.
Both built with DEBUG=TRUE, OPT=TRUE, MPI=TRUE, USE_PETSC=TRUE

```

cd $BISICLES_HOME
git clone git@github.com:ggslc/unicicles.git
cd unicicles

# Build glimmer-cism
cd glimmer-cism
./bootstrap
cd ..
mkdir glimbike-parallel
cd glimbike-parallel

BIKE=$BISICLES_HOME/bisicles-uob/
BIKE_CONFIG=2d.Linux.64.mpiCC.mpif90.DEBUG.OPT.MPI.PETSC
HDF5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/
NETCDF=/usr

PYTHON=python3 FC=gfortran FCFLAGS="-fno-range-check -ffree-line-length-0 -DBISICLES_CDRIVER -DNO_RESCALE -g -I$BIKE/code/src " LDFLAGS="-L$BIKE/code/lib -lBisicles$BIKE_CONFIG -lChomboLibs$BIKE_CONFIG -lpython3.10 -L$HDF5 -lhdf5 -lz " ../glimmer-cism/configure --with-netcdf=$NETCDF --with-hdf5=$HDF5 --prefix=$PWD --disable-python
make
make install

# build wrappers/ukesm-ice_NETCDF
cd wrappers/ukesm-ice_NETCDF
make clean -f Makefile.ubuntu22.04_opt_mpi_petsc
make -f Makefile.ubuntu22.04_opt_mpi_petsc

```




### Ubuntu 20.04 (serial, debug)

```
cd unicicles

# Build glimmer-cism
cd glimmer-cism
./bootstrap
cd ..
mkdir glimbike-serial
cd glimbike-serial
BIKE=$BISICLES_HOME/bisicles-uob/
HDF5=/usr/lib/x86_64-linux-gnu/hdf5/serial/
FC=gfortran FCFLAGS="-fno-range-check -ffree-line-length-0 -DBISICLES_CDRIVER -DNO_RESCALE -g -I$BIKE/code/src " LDFLAGS="-L$BIKE/code/lib -lBisicles2d.Linux.64.g++.gfortran.DEBUG -lChomboLibs2d.Linux.64.g++.gfortran.DEBUG -lpython3.8 -L$HDF5 -lhdf5 -lz " ../glimmer-cism/configure --with-netcdf=/usr --with-hdf5=$HDF5 --prefix=$PWD
make
make install


```

## Testing Notes

### Using TerraFIRMA data

To validate whether pyglint is producing the same integrated SMB over the ice sheets as the Glint code does, we can use the TerraFIRMA data. The SMB as calculated by Glint and passed to BISICLES (for simulations that had active ice sheet components) has been computed using BISICLES filetools.

If we give pyglint the atmos*icecouple.nc files that are produced by the UM and coupling code (the files passed to Glint) we should be able to reproduce the same integrated SMB over the ice sheets as Glint does (allowing for small differences due to different interpolation methods, etc).

The best TerraFIRMA simulation to use for this is labelled u-cs568, a PI-control experiment with interactive ice sheets.

### Testing setup

The code to calculate integracted SMB from the atmos*icecouple.nc files is in the test_TerraFIRMA.py file.

Different versions of pyglint are in different git branches of this repository. To test a particular version of pyglint, checkout the relevant branch and run the test_TerraFIRMA.py file.

The different configuration options are set out in the config/config*.yaml files. To test a particular configuration, pass the relevant config file as an argument to the test_TerraFIRMA.py file. For example:

```bash
> python test_TerraFIRMA.py --config config/cs568_AIS.yaml
```

The outputs from the test_TerraFIRMA.py file are saved in csv files in the outputs/ directory, with a timestamped directory name for each test run. This directory also  contains metadata about the test run, including the git branch and commit of pyglint that was used, and the config file used for the test run.

## Branch ordering

Starting from the 'main' branch, which was the version of pyglint when starting the TerraFIRMA comparisons, the branches build on each other as follows:
- main: original version of pyglint when starting the TerraFIRMA comparisons
- all_lat_lon: removes the use of ilo:ihi etc. which removed the final lat/lon
- topo_vals: changes the topo_mid and topo_max values to match those (we think) used in Glint
- conservation: makes changes to coupling.py:glint_conservation_adjust to implement the total atmosphere column conservation in glint_downscale_gcm
- units: change to unit conversion, more to do to understand whether this is right or not 


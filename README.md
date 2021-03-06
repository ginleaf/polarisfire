# Info:

Collection of scripts to run POLARIS on FIRE simulations for calculating synchrotron emission

## Setup:
Make sure you install polaris https://portia.astrophysik.uni-kiel.de/polaris/ 

Edit your `~/.bashrc` with the following lines at the end to add polaris and polarisfire to your $PYTHONPATH:
*    `export POLARISPATH=/path/to/polaris" `
*    `export POLARISFIREPATH=/path/to/polarisfire" `
*    `export PYTHONPATH=$PYTHONPATH:$POLARISPATH:$POLARISFIREPATH" `

Then run `bash` to refresh and update environment variables.

## Dependencies:

*   `polaris`
*   `h5py`
*   `yt`
*   `unyt`
*   `numpy`
*   `matplotlib`
*   `astropy` if you want to visualize the output fits files of polaris

## Contents:	

1. `mock_example/`: example code to make a grid with constant B field and CR density for testing purposes

* use `make_mock_yt_cube.py` to read in a simulation snapshot and modify its contents to have the desired values of B, ncr etc
* use `converter_octree_fakecubeyt.py` to make the octree

Note: octree converter script is slightly different than that in process_sims since we don't need to recompute n_cr

2. `process_sims/`: scripts to preprocess simulation snapshots and run polaris

* `convert_FIRE_cubes_Suoqing_newyt.py`: holds functions to read Gizmo output and write to grid format. The user can select to zoom within a specified radius and the resulting number of pixels
* `cube_gen.py`: creates a gridded cube (hdf5) using the functions in `convert_FIRE_cubes_Suoqing_newyt.py`
* `cube_gen.sbatch`: batch script to run on wheeler - sends email when done (change email address in script for this)
* `converter_octree_snapshot.py`: converts hdf5 grid file to octree format for input to Polaris        
* `run_polaris.py`: python wrapper that creates command file for Polaris and runs Polaris
* `plot_data.py`: example script to plot fits file outputs of polaris

Run scripts in this order:

1) `python3.7 cube_gen.py` (or the corresponding batch script) to make a gridded cube
2) `python3.7 converter_octree_snapshot.py` to convert that cube to octree
3) `python3.7 run_polaris.py` to run Polaris 

Then you may copy the output fits files (usually in results/.../data/) to local pc
and plot output with `plot_data.py`

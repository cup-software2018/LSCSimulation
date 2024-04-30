# LSCSimulation

## clone and build
	>> git clone https://github.com/cup-software2018/LSCSimulation.git
 	>> cd LSCSimulation
	>> mkdir build; cd build; 
 	>> cmake .. -DCMAKE_INSTALL_PREFIX=[installation directory]
 	>> make -j[NCPU]; make install

## Execution
just in the build directory

	>> ./LSCSim/lscsim
	Usage: LSCSim [-n # of event] [-o output] [-f macro]
              [-g geometry] [-p pmt_position data] [-m material] [-v macro]

where geometry, material, and PMT position data are in LSCSim/data in the source directory and example macro files are in LSCSim/mac.

You can start from the installation directory, but you have to set LD_LIBRARY_PATH before running.

	>> export LD_LIBRARY_PATH=[installation directory/lib]:$LD_LIBRARY_PATH

In the installation directory, you can find mac and data directories.

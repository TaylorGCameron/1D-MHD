# 1D-MHD

This is code that pulls ACE data from CDAWeb, evolves it to the bowshock using the AMRVAC code, then produces a time-series at the bowshock. 

The code runs in Python2.7 and Fortran.

This code requires ai.cdas to pull the data, and an installation of AMRVAC to simulate the data. Because of this, it has to run in linux.

Parameters are set in the amrvac.par file.

To run it, just run the run_model script.

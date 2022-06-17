# lava_2d_py

Hello and thanks for your interest in this project! 

This is the core python code needed to run the lava flow propagation model described in Hyman et al. (2022). 

The main architecture is as follows:

> Most of the model is contained within sim.py
> input.py is for input parameters and pointers to necessary input files such as a DEM (looking for a GeoTIFF), a directory with vent files (all txt), and a directory wher to dump the output file (will be NetCDF) 
> The model can be run from a bash shell with: python input.py 


If something isn't working for you then please let me know.

Cheers!
Dave Hyman

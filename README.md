# lava_2d_py

Hello and thanks for your interest in this project! 
This is the core python code needed to run the lava flow propagation model described in Hyman et al. (2022, in review: JGR - Solid Earth). 
There will likely be a name change to this code, so please consider "lava_2d_py" a placeholder.  

**Please note: This model is not yet peer reviewed, so use at your own risk and discretion**

The main architecture is as follows:

> Most of the model is contained within sim.py

> input.py is for input parameters and pointers to necessary input files such as a DEM (looking for a GeoTIFF), a directory with vent files (all txt), and a directory wher to dump the output file (will be NetCDF) 
> The model can be run from a bash shell with: python input.py 

> rheo.py and therm.py contain submodules for the thermorheologic model, vents.py generates the source term, topo.py ingests and prepares the DEM as well as freezes stagnant, cold lava cells.  

If something isn't working for you then please let me know.

Cheers!
Dave Hyman

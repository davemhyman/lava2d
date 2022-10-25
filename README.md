# Lava2d

Hello and thanks for your interest in this project! 
This is the core python code needed to run the lava flow propagation model described in Hyman et al. (2022, JGR - Solid Earth). 

**Please note: This model is intended for fluid dynamic research purposes only and does not imply endorsement by the U.S. Government. Modeling results do not represent forecasts, official or otherwise.**

The main architecture is as follows:

> Most of the model is contained within sim.py

> input.py is for input parameters and pointers to necessary input files such as a DEM (looking for a GeoTIFF), a directory with vent files (all txt), and a directory wher to dump the output file (will be NetCDF). 

> The model can be run from a bash shell with: python input.py 

> rheo.py and therm.py contain submodules for the thermorheologic model, vents.py generates the source term, topo.py ingests and prepares the DEM as well as freezes stagnant, cold lava cells.  


I have supplied an example vent file in the "example_vents" directory. If something isn't working for you then please let me know.

Cheers,
Dave Hyman


References:

Hyman, D. M. R., Dietterich, H. R., and Patrick, M. R. (2022). Toward next-generation lava flow forecasting: Development of a fast, physics-based lava propagation model. Journal of Geophysical Research: Solid Earth, 127(10):e2022JB024998.


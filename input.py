import sim
#
#-------------------------------------------------------------------------------
sim.set_topo( # set DEM info
    path_to_dem_file    = ('C:\\Users\\...\\dem.tif'),
    Lon_SRC             = 0.00, # source longitude
    Lat_SRC             = 0.00,    # source latitude
    Lon_LowerLeft       = -0.00, # bounding box: lower-left longitude
    Lat_LowerLeft       = -0.00,   # bounding box: lower-left latitude
    Lon_UpperRight      = 0.00, # bounding box: upper-right longitude
    Lat_UpperRight      = 0.00,   # bounding box: upper-right latitude
    fill_level          = 0.0,     # fill topography up to fill_level (m.a.s.l.)
    dx_desired          = 40,          # meters
    smooth_n_times      = 0
    )

#
#-------------------------------------------------------------------------------
sim.set_source( # set vent/fissure info
    path_to_vent_files      = 'C:\\Users\\...\\directory_containing_vent_files\\'
    ) # all vent files must be named vent_01.txt, vent_02.txt, etc

#
#-------------------------------------------------------------------------------
sim.set_vent_props( # set lava properties @ vents
    temperature_vent_C  = 1100, # deg C
    viscosity_vent      = 100, # Pa s
    cryst_vent          = 0.00  # erupted crystal fraction
    )

#
#-------------------------------------------------------------------------------
sim.set_lava_props( # set lava properties throughout
    liquid_density      = 2700,   # kg m-3
    porosity            = 0.4,
    lava_specific_heat  = 1500,   # J kg-1 K-1
    lava_diffusivity    = 5e-7, # m2 s-1
    lava_conductivity   = None,    # W m-1 K-1
    lava_emissivity     = 0.95
    )

#
#-------------------------------------------------------------------------------
sim.set_rheo( # set rheological properties
    yield_strength_crust    = 0, # Pa
    T_core_T_vent_equal     = True  # core temperature equals vent temperature
    )

#
#-------------------------------------------------------------------------------
sim.set_ambient( # set ambient properties
    atm_temperature     = 300, # K
    h_conv              = 50,   # m s-1
    rainfall            = 0,   # m s-1
    ground_temperature  = 300  # K
    )

#
#-------------------------------------------------------------------------------
sim.set_numerics( # set numerical method details
    method = 'OHA',
    efficiency_min      = 10,   # minmum allowable model-clock ratio
    efficiency_max      = 10000, # maximum allowable model-clock ratio
    cfl_max             = 0.5,
    dt_max              = 10.,
    fraction_to_freeze  = 1.0,  # fraction of freezable lava
    tiny_flow           = 0.1,  # min thickness of lava (m)
    )

#
#-------------------------------------------------------------------------------
sim.set_runtime(
    max_iter = None, #one of them can be none, so the default value for both is none, if both none run till killed
    max_time =  12*3600 # s
    )
#
#-------------------------------------------------------------------------------
sim.set_output(
    path_out =  ('C:\\Users\\...\\directory_for_output\\')
    )
#
#-------------------------------------------------------------------------------
#start simulation
sim.run()
#
#-------------------------------------------------------------------------------

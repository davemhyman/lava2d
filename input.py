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
    dx_desired          = 0,          # meters
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
    temperature_vent_C  = 0.0, # deg C
    viscosity_vent      = 0.0, # Pa s
    cryst_vent          = 0.0  # erupted crystal fraction
    )

#
#-------------------------------------------------------------------------------
sim.set_lava_props( # set lava properties throughout
    liquid_density      = 0.0,   # kg m-3
    porosity            = 0.0,
    lava_specific_heat  = 0.0,   # J kg-1 K-1
    lava_diffusivity    = 0.0, # m2 s-1
    lava_conductivity   = None,    # W m-1 K-1
    lava_emissivity     = 0.0
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
    atm_temperature     = 0.0, # K
    h_conv              = 0.0,   # m s-1
    rainfall            = 0.0,   # m s-1
    ground_temperature  = 0.0  # K
    )

#
#-------------------------------------------------------------------------------
sim.set_numerics( # set numerical method details
    method = 'OHA',
    efficiency_min      = 0.0,   # minmum allowable model-clock ratio
    efficiency_max      = 0.0, # maximum allowable model-clock ratio
    cfl_max             = 0.0,
    dt_max              = 0.0,
    fraction_to_freeze  = 0.0,  # fraction of freezable lava
    tiny_flow           = 0.0,  # min thickness of lava (m)
    )

#
#-------------------------------------------------------------------------------
sim.set_runtime(
    max_iter = None, #one of them can be none, so the default value for both is none, if both none run till killed
    max_time =  0.0 # s
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

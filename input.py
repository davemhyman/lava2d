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
sim.set_init( # set initialization type and file
    init_type = None,  # None or 'prior_model' currently supported
    init_file = None  # only specify if init_type = 'prior_model'
    )
#
#-------------------------------------------------------------------------------
sim.set_vent_props( # set lava properties @ vents
    temperature_vent_C  = 1150, # deg C
    viscosity_melt_vent = 100, # Pa s
    cryst_vent          = 0.0  # erupted crystal fraction
    )
#
#-------------------------------------------------------------------------------
sim.set_lava_props( # set lava properties throughout
    liquid_density      = 2700,   # kg m-3
    porosity            = 0.0,
    lava_specific_heat  = 1500,   # J kg-1 K-1
    lava_diffusivity    = 5e-7, # m2 s-1
    lava_conductivity   = None,    # W m-1 K-1
    lava_emissivity     = 0.95
    )
#
#-------------------------------------------------------------------------------
# using Avrami n = 4
sim.set_rheo( # set rheological properties
    phi_max                 = 0.6, # max crystal packing fraction (could be 0.6 e.g., Marsh 1981; Pinkerton and Stevenson 1992), # could be higher (Cashman et al., 1999)
    phi_inf                 = 0.0,
    max_cryst_rate          = 0.0e0, # s-1 # max crystalization rate: max d(phi)/dt
    yield_strength_crust    = 0.0e0, # Pa
    T_core_T_vent_equal     = True  # core temperature equals vent temperature
    )
#
#-------------------------------------------------------------------------------
sim.set_ambient( # set ambient properties
    ground_temperature  = 300, # K
    atm_temperature     = 300, # K
    h_conv              = 50,   # W m-2 K-1
    )
#
#-------------------------------------------------------------------------------
sim.set_numerics( # set numerical method details
    efficiency_min      = 0, # minmum allowable model-clock ratio
    efficiency_max      = 10000,  # maximum allowable model-clock ratio
    cfl_max             = 0.5,
    dt_max              = 5.,
    fraction_to_freeze  = 1.0,  # fraction of freezable lava per time-step
    tiny_flow           = 0.1,  # min thickness of lava (m)
    )
#
#-------------------------------------------------------------------------------
sim.set_runtime(
    max_iter = None, #one of them can be none, so the default value for both is none, if both none run till killed
    max_time_hr = 1*24, # hr
    out_times = [6., 12., 18.], # hr, list of intermediate output times
    run_to_ocean = True
    )
#
#-------------------------------------------------------------------------------
sim.set_output( # where to store out.nc?
    path_out =  ('C:\\Users\\...\\directory_for_output\\')
    )
#
#-------------------------------------------------------------------------------
sim.set_source( # set vent/fissure info: where is vent_nn.txt located?
    path_to_vent_files      = ('C:\\Users\\...\\directory_containing_vent_files\\')
    ) # all vent files must be named vent_01.txt, vent_02.txt, etc
#
#-------------------------------------------------------------------------------
#start simulation
sim.run()
#
#-------------------------------------------------------------------------------
#

import sim
#
#-------------------------------------------------------------------------------
sim.set_topo( # set DEM info
    path_to_dem_file    = ('C:\\Users\\dhyman\\OneDrive - DOI\\Documents\\'
        'CONVERSE\\scenario_SFVF\\SFVolcanicField_DEM\\SFVolcanicField_DEM\\'
        'SFVolcanicFieldDEM_lzw.tif'),
    Lon_SRC             = -111.5507762030, # source longitude
    Lat_SRC             = 35.4221606350,    # source latitude
    Lon_LowerLeft       = -111.64, # bounding box: lower-left longitude
    Lat_LowerLeft       = 35.28,   # bounding box: lower-left latitude
    Lon_UpperRight      = -111.45, # bounding box: upper-right longitude
    Lat_UpperRight      = 35.57,   # bounding box: upper-right latitude
    dx_desired          = 40,          # meters
    smooth_n_times      = 5
    )

#
#-------------------------------------------------------------------------------
sim.set_source( # set vent/fissure info
    path_to_vent_files      = ('C:\\Users\\dhyman\\OneDrive - DOI\\Documents\\'
        'CONVERSE\\scenario_SFVF\\model_results\\syn_eruptive\\'
        'epsiode_02\\run_12\\'),
    discharge_ramp_up       = 600,  # s
    crust_age_at_eruption   = 0     # seconds
    )

#
#-------------------------------------------------------------------------------
sim.set_vent_props( # set lava properties @ vents
    temperature_vent_C  = 1100, # deg C
    viscosity_vent      = 20,  # Pa s
    cryst_vent          = 0.4  # erupted crystal fraction
    )

#
#-------------------------------------------------------------------------------
sim.set_lava_props( # set lava properties throughout
    liquid_density      = 2700,   # kg m-3
    porosity            = 0.1,
    lava_diffusivity    = 5.5e-7, # m2 s-1
    lava_conductivity   = None    # W m-1 K-1
    )

#
#-------------------------------------------------------------------------------
sim.set_rheo( # set rheological properties
    phi_max                 = 0.67, # max crystal packing fraction
    max_cryst_rate          = 0.0e-5, # s-1 # max crystalization rate: max d(phi)/dt
    yield_strength_crust    = 1.0e4, # Pa
    glass_temperature       = None, # K
    T_core_T_vent_equal     = True  # core temperature equals vent temperature
    )

#
#-------------------------------------------------------------------------------
sim.set_ambient( # set ambient properties
    atm_temperature     = 275, # K
    atm_wind            = 10,   # m s-1
    rainfall            = 0,   # m s-1
    ground_temperature  = 275  # K
    )

#
#-------------------------------------------------------------------------------
sim.set_numerics( # set numerical method details
    method = 'OHA',
    efficiency_min      = 20,   # minmum allowable model-clock ratio
    efficiency_max      = 10000, # maximum allowable model-clock ratio
    dt_max              = 10.,
    fraction_to_freeze  = 1.0,  # fraction of freezable lava
    tiny_flow           = 0.1,  # min thickness of lava (m)
    )

#
#-------------------------------------------------------------------------------
sim.set_runtime(
    max_iter = None, #one of them can be none, so the default value for both is none, if both none run till killed
    max_time = 2 * 86400 # s
    )
#
#-------------------------------------------------------------------------------
sim.set_output(
    path_out =  ('C:\\Users\\dhyman\\OneDrive - DOI\\Documents\\'
        'CONVERSE\\scenario_SFVF\\model_results\\syn_eruptive\\'
        'epsiode_02\\run_12\\')
    )
#
#-------------------------------------------------------------------------------
#start simulation
sim.run()
#
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import colorbar as cbar
import matplotlib.colors as c
import time as clock
from geotiff import GeoTiff
import os
from netCDF4 import Dataset

from globals import params as p
from globals import grids as g

import topo
import rheo
import vents
import thermal as therm
import muscl as m
import post



################################################################################
################################################################################
################################################################################
################################################################################
#-------------------------------------------------------------------------------
def restrict_grids():
    # Static:
    g.x = m.restrict(g.x)
    g.y = m.restrict(g.y)
    g.lon = m.restrict(g.lon)
    g.lat = m.restrict(g.lat)
    g.B0 = m.restrict(g.B0)
    p.dx = 2*p.dx
    p.dy = 2*p.dy
    # Dynamic:
    g.h_n = m.restrict(g.h_n)
    g.B_n = m.restrict(g.B_n)
    g.abs_Usurf = m.restrict(g.abs_Usurf)
    g.t_erupted = m.restrict(g.t_erupted)
    g.t_inundation = m.restrict(g.t_inundation)
#
#-------------------------------------------------------------------------------
def interp_grids():
    g.x = m.interp(g.x)
    g.y = m.interp(g.y)
    g.lon = m.interp(g.lon)
    g.lat = m.interp(g.lat)
    g.B0 = m.interp(g.B0)
    p.dx = 0.5*p.dx
    p.dy = 0.5*p.dy
    # Dynamic:
    g.h_n = m.interp(g.h_n)
    g.B_n = m.interp(g.B_n)
    g.abs_Usurf = m.interp(g.abs_Usurf)
    g.t_erupted = m.interp(g.t_erupted)
    g.t_inundation = m.interp(g.t_inundation)
#
#-------------------------------------------------------------------------------
def set_topo(
    path_to_dem_file = '',
    Lon_SRC = 0,
    Lat_SRC = 0,
    Lon_LowerLeft = -.1,
    Lat_LowerLeft = -.1,
    Lon_UpperRight = .1,
    Lat_UpperRight = .1,
    fill_level = None,
    dx_desired = None,
    smooth_n_times = 0
    ):
    #
    topo.topo_from_DEM(
        path_to_dem_file,
        Lon_SRC, Lat_SRC,
        Lon_LowerLeft, Lat_LowerLeft,
        Lon_UpperRight, Lat_UpperRight,
        fill_level = fill_level
        )
    #
    if dx_desired != None:
        G = dx_desired / p.dx
        #
        if G < 1:
            # Go to finer grids (Interpolation)
            n = int(round(np.log2(1./G)))
            for i in range(n):
                g.x = m.interp(g.x)
                g.y = m.interp(g.y)
                g.lon = m.interp(g.lon)
                g.lat = m.interp(g.lat)
                g.B0 = m.interp(g.B0)
                p.dx = 0.5*p.dx
                p.dy = 0.5*p.dy
                #
            #
        #
        elif G > 1:
            # Go to coarser grids (Restriction)
            n = int(round(np.log2(G)))
            for i in range(n):
                g.x = m.restrict(g.x)
                g.y = m.restrict(g.y)
                g.lon = m.restrict(g.lon)
                g.lat = m.restrict(g.lat)
                g.B0 = m.restrict(g.B0)
                p.dx = 2*p.dx
                p.dy = 2*p.dy
                #
            #
        #
    #
    topo.smooth_topo(smooth_n_times)
    print('| --- DEM Preparation Complete --- |')
#
#-------------------------------------------------------------------------------
def set_source(
    path_to_vent_files = ''
    ):
    #
    p.path_to_vent_files = path_to_vent_files
    vents.read_source_data()
#
#-------------------------------------------------------------------------------
def set_vent_props(
    temperature_vent_C = 1140, # deg C
    viscosity_vent = 100, # Pa s
    cryst_vent = 0, # erupted crystal fraction
    ):
    #
    p.T_vent = temperature_vent_C + 273
    p.viscosity_vent = viscosity_vent
    p.cryst_vent = cryst_vent
#
#-------------------------------------------------------------------------------
def set_lava_props(
    liquid_density = 2700, # kg m-3
    porosity = 0.4,
    lava_conductivity = None, # W m-1 K-1
    lava_diffusivity = 5.5e-7, # m2 s-1
    lava_specific_heat = 1000, # J kg-1 K-1
    lava_emissivity = 0.95,
    ):
    #
    p.lava_density = liquid_density * (1 - porosity)
    p.lava_specific_heat = lava_specific_heat
    if lava_conductivity == None:
        p.lava_diffusivity = lava_diffusivity
        p.lava_conductivity = p.lava_diffusivity * p.lava_density * p.lava_specific_heat
    else:
        p.lava_conductivity = lava_conductivity
        p.lava_diffusivity = p.lava_conductivity / (p.lava_density * p.lava_specific_heat)
    #
    p.lava_emissivity = lava_emissivity
#
#-------------------------------------------------------------------------------
def set_rheo(
    phi_max = 0.67, # max crystal packing fraction
    max_cryst_rate = 0, # s-1
    yield_strength_crust = 0, # Pa
    glass_temperature = None, # K
    T_core_T_vent_equal = True
    ):
    #
    p.phi_max = phi_max
    p.max_cryst_rate = max_cryst_rate
    p.yield_strength_crust = yield_strength_crust
    p.glass_temperature = glass_temperature
    if T_core_T_vent_equal:
        p.core_temperature = p.T_vent
    else:
        # whatever cooling occurs near the vent to feed lava flows
        near_vent_speed = 1 # m s-1
        vent_mixing_depth = 2 # m
        q_rad = p.lava_emissivity * p.stefan_boltzmann * p.core_temperature**4
        initial_cooling_rate = q_rad / (p.lava_density * p.lava_specific_heat * vent_mixing_depth)
        time_to_feed_lava_flows = p.dx / near_vent_speed
        p.core_temperature = p.T_vent - initial_cooling_rate * time_to_feed_lava_flows
#
#-------------------------------------------------------------------------------
def set_ambient(
    atm_temperature = 300, # K
    h_conv = 10, # m s-1
    rainfall = 0, # m s-1
    ground_temperature = 300 # K
    ):
    #
    p.atm_temperature = atm_temperature
    p.h_conv = h_conv
    p.rainfall = rainfall
    p.ground_temperature = ground_temperature
#
#-------------------------------------------------------------------------------
def set_numerics(
    method = 'OHA',
    efficiency_min = 0,
    efficiency_max = 10000,
    cfl_max = 0.5,
    dt_max = 10.,
    fraction_to_freeze = 1.0,
    tiny_flow = 0.1
    ):
    #
    p.method = method
    p.efficiency_min = efficiency_min
    p.efficiency_max = efficiency_max
    p.cfl_max = cfl_max
    p.dt_max = dt_max
    p.fraction_to_freeze = fraction_to_freeze
    p.tiny_flow = tiny_flow
#
#-------------------------------------------------------------------------------
def set_runtime(
    max_iter = None,
    max_time_hr = 0 # s
    ):
    #
    time_limited = max_time_hr != None
    iter_limited = max_iter != None
    if not iter_limited and not time_limited:
        # Error: kill immediately
        p.max_time = 0
        p.max_iter = 0
    elif iter_limited and not time_limited:
        # no time limit
        p.max_time = np.inf
        p.max_iter = max_iter
    elif not iter_limited and time_limited:
        # no iter limit
        p.max_time = max_time_hr*3600
        p.max_iter = np.inf
#
#-------------------------------------------------------------------------------
def run():
    print('| --- Initializing Model Domain --- |')
    #
    # ICs: scalars
    vol_erupted = 0.0
    n = 0
    t_n = 0.0
    grid_level = 0
    dts = []
    step_times = []
    n_nodes = []
    #
    # ICs: fields
    g.t_inundation = np.zeros(g.B0.shape, dtype = g.B0.dtype) + t_n
    g.t_erupted = np.zeros(g.B0.shape, dtype = g.B0.dtype)
    g.abs_Usurf = np.zeros(g.B0.shape, dtype = g.B0.dtype)
    g.h_n = np.zeros(g.B0.shape, dtype = g.B0.dtype)
    g.B_n = g.B0.copy()
    #
    buffer = 3 # theoretically, cfl limits buffer to 1
    bounds = m.subset_bounds(t_n, buffer_cells = 2*buffer)
    x, y, h_n, B_n, te_n, ti_n, abs_U_s = m.subset_grids(bounds)
    #I = np.isfinite(x)
    #
    q_n = vents.source_term(x,y,t_n)
    if q_n.max() == 0:
        dt = p.dt_max
    else:
        dt = 2 * p.tiny_flow / q_n[q_n>0].mean()
    #
    print('| --- Beginning Model Evolution --- |')
    t0 = clock.time()
    #
    while (n < p.max_iter) and (t_n < p.max_time):
        #
        t_step_start = clock.time()
        #
        # ----------------------------------------------------------------------
        # Generate Subset
        bounds = m.subset_bounds(t_n, buffer_cells = 2*buffer)
        x, y, h_n, B_n, te_n, ti_n, abs_U_s = m.subset_grids(bounds)
        q_n = vents.source_term(x, y, t_n)
        I = m.active_indicator(h_n, q_n, buffer)
        # ----------------------------------------------------------------------
        # MUSCL step
        #
        # specify separate methods:
        # physics_quality: {'operations', 'research'} (O, R)
        # h_update: {'full', 'parabolic_substeps', 'hyperbolic_approx'} (F, S, H)
        # te_update: {'advection', 'from_inundation'} (A, I)
        #
        #
        if p.method == 'RFA':
            h_next, te_next, abs_U_s, B_next, dt, dt_type, substeps  = m.step_KT2000_FD2_RFA(h_n.copy(), B_n.copy(), te_n.copy(), ti_n.copy(), q_n, abs_U_s, t_n, I)
        elif p.method == 'OSA':
            h_next, te_next, abs_U_s, B_next, dt, dt_type, substeps  = m.step_KT2000_FD2_OSA(h_n.copy(), B_n.copy(), te_n.copy(), ti_n.copy(), q_n, abs_U_s, t_n, I)
        elif p.method == 'OSI':
            h_next, te_next, abs_U_s, B_next, dt, dt_type, substeps  = m.step_KT2000_FD2_OSI(h_n.copy(), B_n.copy(), te_n.copy(), ti_n.copy(), q_n, abs_U_s, t_n, I)
        elif p.method == 'OHA':
            h_next, te_next, abs_U_s, B_next, dt, dt_type, substeps  = m.step_KT2000_FD2_OHA(h_n.copy(), B_n.copy(), te_n.copy(), ti_n.copy(), q_n, abs_U_s, t_n, I)
        elif p.method == 'OHI':
            h_next, te_next, abs_U_s, B_next, dt, dt_type, substeps  = m.step_KT2000_FD2_OHI(h_n.copy(), B_n.copy(), te_n.copy(), ti_n.copy(), q_n, abs_U_s, t_n, I)
        else:
            print('| --- INVALID CHOICE OF METHOD --- |')
            h_next = h_n.copy(); B_next = B_n.copy(); te_next = te_n.copy()
        #
        # ----------------------------------------------------------------------
        # Update Scalars
        vol_erupted += dt * np.sum(q_n) * p.dx * p.dy
        t_n += dt
        n+=1
        #
        # ----------------------------------------------------------------------
        # Re-embed grids:
        r_start, r_stop, c_start, c_stop = bounds
        g.h_n[r_start:r_stop, c_start:c_stop] = h_next # re-embed
        g.B_n[r_start:r_stop, c_start:c_stop] = B_next # re-embed
        g.abs_Usurf[r_start:r_stop, c_start:c_stop] = abs_U_s # re-embed
        g.t_erupted[r_start:r_stop, c_start:c_stop] = te_next # re-embed
        g.t_inundation = m.update_inundation_times(g.t_inundation, g.h_n, t_n, dt) # Backward Searching (t_n already updated)
        #
        # ----------------------------------------------------------------------
        t_step_stop = clock.time()
        step_times.append(t_step_stop - t_step_start)
        n_nodes.append(I.sum())
        dts.append(dt)
        #
        # ----------------------------------------------------------------------
        # every n%10 = 0 (ie: n = 10, 20, 30 ...)
        if n%10 == 0:
            print('|  Step (sub): {} ({})  |  dx: {} m  |  Grid: {} x {}  |  dt (type): {:.03f} s ({})  |  t: {:.02f} hr |'.format(n,substeps, p.dx, h_n.shape[0],h_n.shape[-1], dt, dt_type, t_n/3600.))
        # Decide on Grid Transfer
        #
        if (n >= 100) and (n%100) == 0:
            model_clock_last10 = np.sum(dts[-100:]) / np.sum(step_times[-100:])
            if model_clock_last10 < p.efficiency_min:
                grid_level += 1
                restrict_grids()
            elif model_clock_last10 > p.efficiency_max:
                grid_level -= 1
                interp_grids()
            # else: do nothing
        #
        # ----------------------------------------------------------------------
    #
    # Simulation Ended
    dB = g.B_n - g.B0
    vol_error = np.sum(g.h_n + dB) * p.dx * p.dy / vol_erupted - 1
    #
    t1 = clock.time()
    model_clock_ratio = (t_n-0)/(t1-t0)
    print('|  Step (sub): {} ({})  |  dx: {} m  |  Grid: {} x {}  |  dt (type): {:.03f} s ({})  |  t: {:.02f} hr |'.format(n,substeps, p.dx, h_n.shape[0],h_n.shape[-1], dt, dt_type, t_n/3600.))
    print('|  Clock Time: {:.03f} s |  Model-Clock Ratio : {:.03f}  |  Relative Volume Error : {:.03}  |'.format(t1 - t0, model_clock_ratio, vol_error))
    #
    # Postprocessing:
    print('| --- Postprocessing --- |')
    post.make_all_physics_grids(t_n)
    # Write Output:
    print('| --- Writing Output --- |')
    sim_metadata = [dts, step_times, n_nodes, model_clock_ratio, t1-t0, vol_error, dt, t_n]
    write_nc(sim_metadata)
#
#-------------------------------------------------------------------------------
def set_output(
    path_out = ''
    ):
    #
    # if does not exist, make it
    path_out = os.path.join(path_out, '')
    os.makedirs(path_out, exist_ok = True)
    #
    p.path_out = path_out
#
#-------------------------------------------------------------------------------
def write_nc(sim_metadata):
    #
    dts, step_times, n_nodes, model_clock_ratio, wall_duration, vol_error, dt, t_n = sim_metadata
    dB = g.B_n - g.B0
    h_total = g.h_n + dB
    I = h_total > p.tiny_flow
    advance_rate = np.zeros(g.x.shape)
    abs_grad_ti = rheo.slope_SAN(g.t_inundation, I)
    advance_rate[I] = 1. / (abs_grad_ti + p.pos_eps)
    #
    file = os.path.join(p.path_out, 'out.nc')
    # --------------------------------------------------------------------------
    # Overwrite file
    #
    # --------------------------------------------------------------------------
    if os.path.isfile(file): os.remove(file)
    print('-------------------------------------------------')
    print(clock.ctime(clock.time()))
    print('GENERATING NETCDF FILE:')
    print(file)
    #
    # --------------------------------------------------------------------------
    # OPEN NEW DATASET
    dset = Dataset(file, 'w')
    #
    g_meta      = dset.createGroup('METADATA')
    g_data      = dset.createGroup('DATA')
    g_data_par  = g_data.createGroup('PARAMS')
    g_data_phy  = g_data.createGroup('PHYSICS')
    g_data_haz  = g_data.createGroup('HAZARDS')
    #
    # DIMENSIONS
    ny,nx = g.B0.shape
    n_x         = g_data.createDimension('cols',nx)
    n_y         = g_data.createDimension('rows',ny)
    n_i         = g_meta.createDimension('iters', len(dts))
    scalar      = g_meta.createDimension('parameter', None)
    scalar      = g_data_par.createDimension('parameter', None)
    #
    # VARIABLES
    dset_dts        = g_meta.createVariable('all_timesteps', np.float32, ('iters',))
    dset_clock      = g_meta.createVariable('all_clock', np.float32, ('iters',))
    dset_cells      = g_meta.createVariable('all_num_cells', np.float32, ('iters',))
    dset_tmax       = g_meta.createVariable('model_duration','f4', ('parameter',))
    dset_wall       = g_meta.createVariable('wall_runtime','f4', ('parameter',))
    dset_vol_err    = g_meta.createVariable('volume_error','f4', ('parameter',))
    #
    dset_x          = g_data.createVariable('x', np.float32, ('rows','cols')) # (ny,nx)
    dset_y          = g_data.createVariable('y', np.float32, ('rows','cols')) # (ny,nx)
    dset_lon        = g_data.createVariable('lon', np.float32, ('rows','cols')) # (ny,nx)
    dset_lat        = g_data.createVariable('lat', np.float32, ('rows','cols')) # (ny,nx)
    dset_B_init     = g_data.createVariable('topo_init', np.float32, ('rows','cols')) # (ny,nx)
    #
    dset_T_vent                 = g_data_par.createVariable('vent_temperature','f4', ('parameter',))
    dset_mu_vent                = g_data_par.createVariable('vent_viscosity','f4', ('parameter',))
    dset_cryst_vent             = g_data_par.createVariable('cryst_vent','f4', ('parameter',))
    dset_lava_density           = g_data_par.createVariable('lava_density','f4', ('parameter',))
    dset_lava_specific_heat     = g_data_par.createVariable('lava_specific_heat','f4', ('parameter',))
    dset_lava_conductivity      = g_data_par.createVariable('lava_conductivity','f4', ('parameter',))
    dset_lava_diffusivity       = g_data_par.createVariable('lava_diffusivity','f4', ('parameter',))
    dset_lava_emissivity        = g_data_par.createVariable('lava_emissivity','f4', ('parameter',))
    dset_phi_max                = g_data_par.createVariable('phi_max','f4', ('parameter',))
    dset_max_cryst_rate         = g_data_par.createVariable('max_cryst_rate','f4', ('parameter',))
    dset_yield_strength_crust   = g_data_par.createVariable('yield_strength_crust','f4', ('parameter',))
    dset_core_temperature       = g_data_par.createVariable('core_temperature','f4', ('parameter',))
    dset_atm_temperature        = g_data_par.createVariable('atm_temperature','f4', ('parameter',))
    dset_h_conv                 = g_data_par.createVariable('h_conv','f4', ('parameter',))
    dset_rainfall               = g_data_par.createVariable('rainfall','f4', ('parameter',))
    dset_ground_temperature     = g_data_par.createVariable('ground_temperature','f4', ('parameter',))
    dset_efficiency_min         = g_data_par.createVariable('efficiency_min','f4', ('parameter',))
    dset_efficiency_max         = g_data_par.createVariable('efficiency_max','f4', ('parameter',))
    dset_cfl_max                = g_data_par.createVariable('cfl_max','f4', ('parameter',))
    dset_dt_max                 = g_data_par.createVariable('dt_max','f4', ('parameter',))
    dset_fraction_to_freeze     = g_data_par.createVariable('fraction_to_freeze','f4', ('parameter',))
    dset_tiny_flow              = g_data_par.createVariable('tiny_flow','f4', ('parameter',))
    #
    dset_extent     = g_data_haz.createVariable('extent', np.float32, ('rows','cols')) # (ny,nx)
    dset_arrival    = g_data_haz.createVariable('arrival_time', np.float32, ('rows','cols')) # (ny,nx)
    dset_advance    = g_data_haz.createVariable('advance_rate', np.float32, ('rows','cols')) # (ny,nx)
    #
    dset_B          = g_data_phy.createVariable('topo_dyn', np.float32, ('rows','cols')) # (ny,nx)
    dset_h          = g_data_phy.createVariable('lava_thickness_dyn', np.float32, ('rows','cols')) # (ny,nx)
    dset_h_tot      = g_data_phy.createVariable('lava_thickness_total', np.float32, ('rows','cols')) # (ny,nx)
    dset_dB         = g_data_phy.createVariable('basal_change', np.float32, ('rows','cols')) # (ny,nx)
    dset_t_ti       = g_data_phy.createVariable('time_since_inundation', np.float32, ('rows','cols')) # (ny,nx)
    dset_t_te       = g_data_phy.createVariable('time_since_erupted', np.float32, ('rows','cols')) # (ny,nx)
    dset_abs_Usurf  = g_data_phy.createVariable('surface_speed', np.float32, ('rows','cols')) # (ny,nx)
    dset_Ts         = g_data_phy.createVariable('surface_temperature', np.float32, ('rows','cols')) # (ny,nx)
    #
    # VARIABLE ATTRIBUTES
    dset_dts.description        = 'timestep size for each iteration'
    dset_clock.description      = 'wall clock time for each iteration'
    dset_cells.description      = 'number of active cells for each iteration'
    dset_tmax.description       = 'model duration'
    dset_wall.description       = 'wall clock duration'
    dset_vol_err.description    = 'lava volume relative error: deposit vs. erupted'
    #
    dset_x.description      = 'distance east of source'
    dset_y.description      = 'distance north of source'
    dset_lon.description    = 'longitude'
    dset_lat.description    = 'latitude'
    dset_B_init.description = 'initial DEM'
    #
    dset_T_vent.description                 = 'vent temperature (well-mixed assumption)'
    dset_mu_vent.description                = 'vent viscosity (liquid)'
    dset_cryst_vent.description             = 'erupted crystal fraction'
    dset_lava_density.description           = 'bulk density'
    dset_lava_specific_heat.description     = 'lumped specific heat'
    dset_lava_conductivity.description      = 'lumped thermal conductivity'
    dset_lava_diffusivity.description       = 'lumped thermal diffusivity'
    dset_lava_emissivity.description        = 'lava surface emissivity'
    dset_phi_max.description                = 'max packing crystal fraction'
    dset_max_cryst_rate.description         = 'largest crystallization rate'
    dset_yield_strength_crust.description   = 'yield strength of crustal material'
    dset_core_temperature.description       = 'flow core temperature (insulated)'
    dset_atm_temperature.description        = 'ambient temperature above flow'
    dset_h_conv.description                 = 'total heat transfer coefficient (natural and convective)'
    dset_rainfall.description               = 'rainfall rate'
    dset_ground_temperature.description     = 'ambient temperature below flow'
    dset_efficiency_min.description         = 'min allowable model-clock ratio before coarser grid used'
    dset_efficiency_max.description         = 'max allowable model-clock ratio before finer grid used'
    dset_cfl_max.description                = 'max CFL number (scheme-dependent, must be less than 1)'
    dset_dt_max.description                 = 'max allowable time step'
    dset_fraction_to_freeze.description     = 'fraction of lava thickness to freeze per time step for eligible cells'
    dset_tiny_flow.description              = 'minimum real lava thickness'
    #
    dset_extent.description     = 'lava inundation flag'
    dset_arrival.description    = 'lava arrival time at each cell'
    dset_advance.description    = 'lava advance rate at each cell'
    #
    dset_B.description          = 'dynamic topography = initial DEM + basal change'
    dset_h.description          = 'dynamic lava thickness'
    dset_h_tot.description      = 'total deposit thickness = dynamic lava thickness + basal change'
    dset_dB.description         = 'vertical topographic changes (freezing / remelting)'
    dset_t_ti.description       = 'time since each cell was inundated'
    dset_t_te.description       = 'lava surface travel time from vent to cell'
    dset_abs_Usurf.description  = 'lava surface speed (magnitude)'
    dset_Ts.description         = 'lava surface temperature'

    dset_dts.units      = 's'
    dset_clock.units    = 's'
    dset_cells.units    = '-'
    dset_tmax.units     = 's'
    dset_wall.units     = 's'
    dset_vol_err.units  = '-'
    #
    dset_x.units        = 'm'
    dset_y.units        = 'm'
    dset_lon.units      = 'decimal degrees east'
    dset_lat.units      = 'decimal degrees north'
    dset_B_init.units   = 'm asl'
    #
    dset_T_vent.units               = 'K'
    dset_mu_vent.units              = 'Pa s'
    dset_cryst_vent.units           = '-'
    dset_lava_density.units         = 'kg m-3'
    dset_lava_specific_heat.units   = 'J kg-1 K-1'
    dset_lava_conductivity.units    = 'W m-1 K-1'
    dset_lava_diffusivity.units     = 'm2 s-1'
    dset_lava_emissivity.units      = '-'
    dset_phi_max.units              = '-'
    dset_max_cryst_rate.units       = 's-1'
    dset_yield_strength_crust.units = 'Pa'
    dset_core_temperature.units     = 'K'
    dset_atm_temperature.units      = 'K'
    dset_h_conv.units               = 'W m-2 K-1'
    dset_rainfall.units             = 'm s-1'
    dset_ground_temperature.units   = 'K'
    dset_efficiency_min.units       = '-'
    dset_efficiency_max.units       = '-'
    dset_cfl_max.units              = '-'
    dset_dt_max.units               = 's'
    dset_fraction_to_freeze.units   = '-'
    dset_tiny_flow.units            = 'm'
    #
    dset_extent.units   = '-'
    dset_arrival.units  = 's'
    dset_advance.units  = 'm s-1'
    #
    dset_B.units            = 'm asl'
    dset_h.units            = 'm'
    dset_h_tot.units        = 'm'
    dset_dB.units           = 'm'
    dset_t_ti.units         = 's'
    dset_t_te.units         = 's'
    dset_abs_Usurf.units    = 'm s-1'
    dset_Ts.units           = 'C'
    #
    # ADD VALUES TO VARIABLES
    dset_dts[:]         = np.float32(dts)
    dset_clock[:]       = np.float32(step_times)
    dset_cells[:]       = np.float32(n_nodes)
    dset_tmax[0]        = np.float32(t_n)
    dset_wall[0]        = np.float32(wall_duration)
    dset_vol_err[0]     = np.float32(vol_error)
    #
    dset_x[:]       = np.float32(g.x)
    dset_y[:]       = np.float32(g.y)
    dset_lon[:]     = np.float32(g.lon)
    dset_lat[:]     = np.float32(g.lat)
    dset_B_init[:]  = np.float32(g.B0)
    #
    dset_T_vent[0]                 = np.float32(p.T_vent)
    dset_mu_vent[0]                = np.float32(p.viscosity_vent)
    dset_cryst_vent[0]             = np.float32(p.cryst_vent)
    dset_lava_density[0]           = np.float32(p.lava_density)
    dset_lava_specific_heat[0]     = np.float32(p.lava_specific_heat)
    dset_lava_conductivity[0]      = np.float32(p.lava_conductivity)
    dset_lava_diffusivity[0]       = np.float32(p.lava_diffusivity)
    dset_lava_emissivity[0]        = np.float32(p.lava_emissivity)
    dset_phi_max[0]                = np.float32(p.phi_max)
    dset_max_cryst_rate[0]         = np.float32(p.max_cryst_rate)
    dset_yield_strength_crust[0]   = np.float32(p.yield_strength_crust)
    dset_core_temperature[0]       = np.float32(p.core_temperature)
    dset_atm_temperature[0]        = np.float32(p.atm_temperature)
    dset_h_conv[0]                 = np.float32(p.h_conv)
    dset_rainfall[0]               = np.float32(p.rainfall)
    dset_ground_temperature[0]     = np.float32(p.ground_temperature)
    dset_efficiency_min[0]         = np.float32(p.efficiency_min)
    dset_efficiency_max[0]         = np.float32(p.efficiency_max)
    dset_cfl_max[0]                = np.float32(p.cfl_max)
    dset_dt_max[0]                 = np.float32(p.dt_max)
    dset_fraction_to_freeze[0]     = np.float32(p.fraction_to_freeze)
    dset_tiny_flow[0]              = np.float32(p.tiny_flow)
    #
    dset_extent[:]      = np.float32(I)
    dset_arrival[:]     = np.float32(g.t_inundation)
    dset_advance[:]     = np.float32(advance_rate)
    #
    dset_B[:]           = np.float32(g.B_n)
    dset_h[:]           = np.float32(g.h_n)
    dset_h_tot[:]       = np.float32(h_total)
    dset_dB[:]          = np.float32(dB)
    dset_t_ti[:]        = np.float32(t_n - g.t_inundation)
    dset_t_te[:]        = np.float32(t_n - g.t_erupted)
    dset_abs_Usurf[:]   = np.float32(g.abs_Usurf)
    dset_Ts[:]          = np.float32(g.surface_T_n)
    #
    # WRITE FILE
    dset.close()











#

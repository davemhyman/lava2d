import numpy as np

class params():
    # Stefan-Boltzmann
    stefan_boltzmann = 5.67e-8 # W m-2 K-4
    # Gravity
    g = 9.81 # m s-2
    # Water
    water_density = 1000 # kg m-3
    water_specific_heat = 4190 # J kg-1 K-1
    water_latent_heat_vap = 2.26e6 # J kg-1
    water_boiling_temperature = 373 # K
    #
    pos_eps = 1e-16
    lambda_max = 0.5
    safety_factor = 0.9
    #
    vent_param_splines = []



class grids():
    #
    # Static
    x = np.array([], dtype = np.float64)
    y = np.array([], dtype = np.float64)
    B0 = np.array([], dtype = np.float64)
    vent_grid_mask = np.array([], dtype = np.float64)
    vent_spatial_density = np.array([], dtype = np.float64)
    #
    # Dynamic
    B_n = np.array([], dtype = np.float64)
    h_n = np.array([], dtype = np.float64)
    R_n = np.array([], dtype = np.float64)
    Rs_n = np.array([], dtype = np.float64)
    #
    t_inundation = np.array([], dtype = np.float64)
    t_erupted = np.array([], dtype = np.float64)
    #













#

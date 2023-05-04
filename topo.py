import numpy as np
import rasterio
from rasterio.warp import transform
from rasterio.transform import Affine

from globals import params as p
from globals import grids as g
import muscl as m
import rheo
import thermal as therm

################################################################################
################################################################################
################################################################################
################################################################################

#
############################################################################
#
def topo_from_DEM(
    file,
    Lon_SRC, Lat_SRC,
    Lon_LowerLeft, Lat_LowerLeft,
    Lon_UpperRight, Lat_UpperRight,
    fill_level = None
    ):
    #
    print('| --- Reading DEM Dataset --- |')
    dataset = rasterio.open(file,'r')
    rows,cols = dataset.height, dataset.width
    p.dx = (dataset.bounds.right - dataset.bounds.left) / cols
    p.dy = (dataset.bounds.top - dataset.bounds.bottom) / rows
    #
    # Note: for theoretical topography: lon == x, lat == y
    lons = [Lon_SRC, Lon_LowerLeft, Lon_UpperRight]
    lats = [Lat_SRC, Lat_LowerLeft, Lat_UpperRight]
    xs, ys = transform({'init': 'EPSG:4326'}, dataset.crs, lons, lats)
    x_src, x_LL, x_UR = xs
    y_src, y_LL, y_UR = ys
    #
    width = x_UR - x_LL
    height = y_UR - y_LL
    x_frac = (x_src - x_LL) / width
    y_frac = (y_src - y_LL) / height
    #
    x_levels = int(np.round(np.log2(width/p.dx)))
    y_levels = int(np.round(np.log2(height/p.dy)))
    #
    xvec = dataset.bounds.left + np.arange(cols)*p.dx
    yvec = dataset.bounds.top - np.arange(rows)*p.dy
    x, y = np.meshgrid(xvec,yvec)
    #
    x_source = xvec[np.argmin(abs(xvec - x_src))]
    y_source = yvec[np.argmin(abs(yvec - y_src))]
    #
    nx = 2**x_levels + 1
    ny = 2**y_levels + 1
    #
    nx_off = int(x_frac*nx)
    ny_off = int(y_frac*ny)
    #
    xmin = x_source - nx_off*p.dx
    xmax = x_source + (nx-nx_off)*p.dx
    ymin = y_source - ny_off*p.dy
    ymax = y_source + (ny-ny_off)*p.dy
    #
    c_left = np.argmin(abs(xvec-xmin))
    c_right = np.argmin(abs(xvec-xmax))
    r_top = np.argmin(abs(yvec-ymax))
    r_bottom = np.argmin(abs(yvec-ymin))
    #
    print('| --- Building Modelling Subset --- |')
    dem = dataset.read(1)[r_top:r_bottom,c_left:c_right]
    x_sub = x[r_top:r_bottom,c_left:c_right]
    y_sub = y[r_top:r_bottom,c_left:c_right]
    #
    print('| --- Building Coordinate Transform --- |')
    # need to get lat lon grids:
    lon_sub, lat_sub = transform(dataset.crs, {'init': 'EPSG:4326'}, x_sub.reshape(-1), y_sub.reshape(-1))
    g.lon = np.array(lon_sub).reshape(x_sub.shape)
    g.lat = np.array(lat_sub).reshape(y_sub.shape)
    #
    g.x = (x_sub - x_source).astype(g.x.dtype)
    g.y = (y_sub - y_source).astype(g.x.dtype)
    #
    print('| --- Preparing DEM --- |')
    if fill_level != None:
        ocean = dem <= fill_level
        dem[ocean] = fill_level
        p.dem_fill_level = fill_level
    g.B0 = dem.copy().astype(g.x.dtype)


def smooth_topo(n_times):
    for i in range(n_times):
        print('| --- Smoothing Topography --- |')
        g.B0 = m.local_mean(g.B0)



def freeze(h, B_n, dBdt, ti, phi_S, cryst_core, jeffreys_efficiency, t_n):
    #
    # Implicitly global:
    # phi_S, cryst_core, newtonian_efficiency, t_n
    #
    # global h, B_n, dBdt, ti
    #
    # Freezout:
    # if cell to cell velocity is zero and core is locked or crust is full thickness
    # Locked core and fully crust only casues cell-to-cell zero veolicty when these conditions are true in both cells.
    # Only true at cell ij when also true for E, W, N, S
    #
    T_stiff = rheo.T_stiffened(p.core_temperature) # temp at which visc = 2*visc_core
    eta_inf = therm.erfinv((T_stiff-p.atm_temperature)/(p.core_temperature-p.atm_temperature))
    phi_true_crust = phi_S * (0.88/eta_inf) # 0.88 from Hon et al., 1998
    #
    fully_crust = phi_true_crust >= 1.0 # all of flow thickness is crust
    core_locked = cryst_core >= p.phi_max # if all neighborhood cells true
    stalled = m.local_and_EWNS(jeffreys_efficiency == 0) # EWNS neighbors all stopped
    freezeout = np.logical_and(stalled, np.logical_or(core_locked, fully_crust))
    #
    h_to_freeze = p.fraction_to_freeze * h[freezeout]
    B_n[freezeout] += h_to_freeze
    h[freezeout] -= h_to_freeze
    dBdt[freezeout] = 0
    #







#

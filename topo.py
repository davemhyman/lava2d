import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import colorbar as cbar
import matplotlib.colors as c
#from geotiff import GeoTiff
import rasterio
from rasterio.warp import transform
from scipy import interpolate

from globals import params as p
from globals import grids as g


import muscl as m
import post









################################################################################
################################################################################
################################################################################
################################################################################

#
############################################################################
#
'''
def topo_from_DEM(file, lon_src, lat_src, width, height, x_frac = 0.1, y_frac = 0.1, rot_angle = 0):
    gt = GeoTiff(file)
    rows,cols = gt.tifShape
    crs = gt.crs_code
    x_src, y_src = gt._convert_from_wgs_84(crs, [lon_src,lat_src])
    #
    ((x_left,y_top),(x_right,y_bottom))  = gt.tif_bBox
    p.dx = (x_right - x_left) / cols
    p.dy = (y_top - y_bottom) / rows
    #
    x_levels = int(np.round(np.log2(width/p.dx)))
    y_levels = int(np.round(np.log2(height/p.dy)))
    nx = 2**x_levels + 1
    ny = 2**y_levels + 1
    n_buffer = int((nx**2 + ny**2)**0.5)
    #
    xvec = x_left + np.arange(cols) * p.dx
    yvec = y_top - np.arange(rows) * p.dy
    xs_vec = xvec - x_src
    ys_vec = yvec - y_src
    idx_x_src = np.argmin(abs(xs_vec))
    idx_y_src = np.argmin(abs(ys_vec))
    xs_vec_sub = xs_vec[idx_x_src-n_buffer:idx_x_src+n_buffer]
    ys_vec_sub = ys_vec[idx_y_src-n_buffer:idx_y_src+n_buffer]
    xs_sub, ys_sub = np.meshgrid(xs_vec_sub,ys_vec_sub)
    #
    dem = gt.read()[idx_y_src-n_buffer:idx_y_src+n_buffer, idx_x_src-n_buffer:idx_x_src+n_buffer]
    ocean = dem <= 0.0
    dem[ocean] = 0.0
    #
    # Make output grids (g.x, g.y):
    nx_off = int(x_frac*nx)
    ny_off = int(y_frac*ny)
    xmin = -nx_off*p.dx
    xmax = (nx-nx_off)*p.dx
    ymin = -ny_off*p.dy
    ymax = (ny-ny_off)*p.dy
    c_left = np.argmin(abs(xs_vec_sub-xmin))
    c_right = np.argmin(abs(xs_vec_sub-xmax))
    r_top = np.argmin(abs(ys_vec_sub-ymax))
    r_bottom = np.argmin(abs(ys_vec_sub-ymin))
    g.x = xs_sub[r_top:r_bottom,c_left:c_right].astype(g.x.dtype)
    g.y = ys_sub[r_top:r_bottom,c_left:c_right].astype(g.x.dtype)
    #
    # Angle anti-clockwise from East (mathematical polar angle)
    ROT = rotation_matrix(rot_angle) # 2 x 2 matrix
    rows, cols = dem.shape
    x_all = xs_sub.reshape(-1)
    y_all = ys_sub.reshape(-1)
    position_all = np.concatenate((x_all[None,:], y_all[None,:]), axis = 0) # 2 x N array of column vectors
    position_all_rot = ROT.dot(position_all) # 2 x N array of column vectors
    x_new_all = position_all_rot[0,:]
    y_new_all = position_all_rot[1,:]
    valid = (x_new_all>=xmin-5*p.dx)*(x_new_all<=xmax+5*p.dx)*(y_new_all>=ymin-5*p.dy)*(y_new_all<=ymax+5*p.dy)
    dem_new = interpolate.griddata((x_new_all[valid], y_new_all[valid]), dem.reshape(-1)[valid], (g.x, g.y), method='linear')
    dem_new_corrected = interpolate.griddata((g.x[np.isfinite(dem_new)], g.y[np.isfinite(dem_new)]), dem_new[np.isfinite(dem_new)], (g.x, g.y), method='nearest')
    #
    g.B0 = dem_new_corrected.copy().astype(g.x.dtype)

'''
'''
def topo_from_DEM(file, lon_src, lat_src, width, height, x_frac = 0.5, y_frac = 0.5):
    gt = GeoTiff(file)
    rows,cols = gt.tifShape
    crs = gt.crs_code
    x_src, y_src = gt._convert_from_wgs_84(crs, [lon_src,lat_src])
    #
    ((x_left,y_top),(x_right,y_bottom))  = gt.tif_bBox
    p.dx = (x_right - x_left) / cols
    p.dy = (y_top - y_bottom) / rows
    #
    x_levels = int(np.round(np.log2(width/p.dx)))
    y_levels = int(np.round(np.log2(height/p.dy)))
    #
    xvec = x_left + np.arange(cols) * p.dx
    yvec = y_top - np.arange(rows) * p.dy
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
    dem = gt.read()[r_top:r_bottom,c_left:c_right]
    g.x = (x[r_top:r_bottom,c_left:c_right] - x_source).astype(g.x.dtype)
    g.y = (y[r_top:r_bottom,c_left:c_right] - y_source).astype(g.x.dtype)
    #
    ocean = dem <= 0.0
    dem[ocean] = 0.0
    g.B0 = dem.copy().astype(g.x.dtype)
'''
'''
def topo_from_DEM(file, Lon_SRC, Lat_SRC, Lon_LowerLeft, Lat_LowerLeft, Lon_UpperRight, Lat_UpperRight):
    #
    #gt = GeoTiff(file)
    #rows,cols = gt.tif_shape
    #crs = gt.crs_code
    #((x_left,y_top),(x_right,y_bottom))  = gt.tif_bBox
    #p.dx = (x_right - x_left) / cols
    #p.dy = (y_top - y_bottom) / rows
    #
    #area_box = (Lon_LowerLeft, Lat_UpperRight), (Lon_UpperRight, Lat_LowerLeft)
    #A = gt.read_box(area_box, outer_points=2)
    #lon, lat = gt.get_coord_arrays(area_box, outer_points=2)
    #
    dataset = rasterio.open(file,'r')
    rows,cols = dataset.height, dataset.width
    p.dx = (dataset.bounds.right - dataset.bounds.left) / cols
    p.dy = (dataset.bounds.top - dataset.bounds.bottom) / rows
    #
    x_src, y_src = rasterio.warp.transform({'init': 'EPSG:4326'}, dataset.crs, [Lon_SRC] , [Lat_SRC])
    x_LL, y_LL = rasterio.warp.transform({'init': 'EPSG:4326'}, dataset.crs, [Lon_LowerLeft], [Lat_LowerLeft])
    x_UR, y_UR = rasterio.warp.transform({'init': 'EPSG:4326'}, dataset.crs, [Lon_UpperRight], [Lat_UpperRight])
    #
    width = x_UR - x_LL
    height = y_UR - y_LL
    x_frac = (x_src - x_LL) / width
    y_frac = (y_src - y_LL) / height
    #
    x_levels = int(np.round(np.log2(width/p.dx)))
    y_levels = int(np.round(np.log2(height/p.dy)))
    #
    xvec = x_left + np.arange(cols) * p.dx
    yvec = y_top - np.arange(rows) * p.dy
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
    dem = gt.read()[r_top:r_bottom,c_left:c_right]
    g.x = (x[r_top:r_bottom,c_left:c_right] - x_source).astype(g.x.dtype)
    g.y = (y[r_top:r_bottom,c_left:c_right] - y_source).astype(g.x.dtype)
    #
    ocean = dem <= 0.0
    dem[ocean] = 0.0
    g.B0 = dem.copy().astype(g.x.dtype)
'''


def topo_from_DEM(file, Lon_SRC, Lat_SRC, Lon_LowerLeft, Lat_LowerLeft, Lon_UpperRight, Lat_UpperRight):
    #
    print('| --- Reading DEM Dataset --- |')
    dataset = rasterio.open(file,'r')
    rows,cols = dataset.height, dataset.width
    p.dx = (dataset.bounds.right - dataset.bounds.left) / cols
    p.dy = (dataset.bounds.top - dataset.bounds.bottom) / rows
    #
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
    ocean = dem <= 0.0
    dem[ocean] = 0.0
    g.B0 = dem.copy().astype(g.x.dtype)





def smooth_topo(n_times):
    for i in range(n_times):
        print('| --- Smoothing Topography --- |')
        g.B0 = m.local_mean(g.B0)


def topo_from_theory(slope_ang, dx = 10, dy = 10, x_frac = 0.1, y_frac = 0.1, x_levels = 10, y_levels = 10):
    #
    p.dx = dx
    p.dy = dy
    #
    nx = 2**x_levels + 1
    ny = 2**y_levels + 1
    #
    nx_off = int(x_frac*nx)
    ny_off = int(y_frac*ny)
    #
    xmin = -nx_off*p.dx
    xmax = (nx-nx_off)*p.dx
    ymin = -ny_off*p.dy
    ymax = (ny-ny_off)*p.dy
    #
    dBdx = np.tan(np.radians(slope_ang)).astype(g.x.dtype)
    #
    xvec = (xmin + np.arange(nx) * p.dx).astype(g.x.dtype)
    yvec = (ymax - np.arange(ny) * p.dy).astype(g.x.dtype)
    #
    g.x, g.y = np.meshgrid(xvec, yvec)
    slope = - dBdx * g.x
    g.B0 = slope - slope.min()


def rotation_matrix(ang, units = 'rad'):
    # --------------------------------------------------------------------------
    # Compute 2x2 rotation matrix for a given angle
    #
    # --------------------------------------------------------------------------
    if units == 'rad':
        return np.array([[np.cos(ang), -np.sin(ang)],[np.sin(ang), np.cos(ang)]])
    elif units == 'deg':
        ang_rad = ang * np.pi / 180.
        return np.array([[np.cos(ang_rad), -np.sin(ang_rad)],[np.sin(ang_rad), np.cos(ang_rad)]])
    else:
        print('Specify input angle units: "rad" (default) or "deg".')



def slope_angle(B, dx, dy):
    B_ij = B[1:-1,1:-1] # right
    dBdx_R = (B[1:-1,2:] - B_ij) / dx # x-deriv on right
    dBdx_L = (B_ij - B[1:-1,0:-2]) / dx # x-deriv on left
    dBdy_U = (B[0:-2,1:-1] - B_ij) / dy # y-deriv on upper
    dBdy_D = (B_ij - B[2:,1:-1]) / dy # y-deriv on lower (down)
    slope_RU = np.arctan(np.sqrt(dBdx_R**2 + dBdy_U**2))
    slope_RD = np.arctan(np.sqrt(dBdx_R**2 + dBdy_D**2))
    slope_LU = np.arctan(np.sqrt(dBdx_L**2 + dBdy_U**2))
    slope_LD = np.arctan(np.sqrt(dBdx_L**2 + dBdy_D**2))
    B_slope = np.zeros(B.shape)
    B_slope[1:-1,1:-1] = (slope_RU + slope_RD + slope_LU + slope_LD) / 4. * 180 / np.pi
    return B_slope





def show_slope(B, dz = 10):
    #
    B_contours = np.arange(np.floor(B.min()/dz), np.ceil(B.max()/dz) + 2) * dz
    B_slope = slope_angle(B, p.dx, p.dy)
    #
    fig, ax1 = plt.subplots(figsize=(16,8))
    plt.contourf(g.x, g.y, B_slope, np.arange(0,20.1,0.5), cmap= 'jet');plt.colorbar()
    plt.axis('equal')
    plt.plot(0, 0, 'r^', ms = 10)
    plt.contour(g.x, g.y, B, B_contours, colors = 'grey');plt.show()




def show_topo(B, dz = 10):
    #
    B_contours = np.arange(np.floor(B.min()/dz), np.ceil(B.max()/dz) + 2) * dz
    new_cmap = post.truncate_colormap(plt.get_cmap('terrain'), 0.25, 0.9)
    #
    fig, ax1 = plt.subplots(figsize=(16,8))
    plt.contourf(g.x, g.y, B, B_contours, cmap= new_cmap)
    plt.plot(0, 0, 'r^', ms = 10)
    plt.axis('equal')
    cbar();plt.show()




def freeze():
    #
    # Implicitly global:
    # phi_S, cryst_core, newtonian_efficiency, t_n
    #
    global h, B_n, dBdt, ti
    #
    # Freezout:
    # if cell to cell velocity is zero and core is locked or crust is full thickness
    # Locked core and fully crust only casues cell-to-cell zero veolicty when these conditions are true in both cells.
    # Only true at cell ij when also true for E, W, N, S
    #
    fully_crust = phi_S >= 1.0 # if all neighborhood cells true
    core_locked = cryst_core >= p.phi_max # if all neighborhood cells true
    stalled = m.local_and_EWNS(newtonian_efficiency == 0) # EWNS neighbors all stopped
    freezeout = np.logical_and(stalled, np.logical_or(core_locked, fully_crust))
    #
    h_to_freeze = p.fraction_to_freeze * h[freezeout]
    B_n[freezeout] += h_to_freeze
    h[freezeout] -= h_to_freeze
    dBdt[freezeout] = 0
    ti[freezeout] = t_n
    # R, Rs, abs_U_s all zero already in freezout region
    # no return needed







#

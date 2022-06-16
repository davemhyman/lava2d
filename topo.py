import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import colorbar as cbar
import matplotlib.colors as c
#from geotiff import GeoTiff
import rasterio
from rasterio.warp import transform
from rasterio.transform import Affine
from scipy import interpolate

from globals import params as p
from globals import grids as g


import muscl as m
import rheo
import thermal as therm
import post

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
    g.B0 = dem.copy().astype(g.x.dtype)


def make_inclined_plane(dip, azimuth, x_UpperLeft, y_UpperLeft, nx, ny, dxy, path_to_dir):
    # Make inclined plane dem
    # azimuth: direction down (degrees E of N)
    # dip: how steep (degrees down from horizontal)
    # lat,lon will be equal to x,y
    #
    x = x_UpperLeft + np.arange(nx)*dxy
    y = y_UpperLeft - np.arange(ny)*dxy
    X,Y = np.meshgrid(x,y)
    #
    azir = np.radians(azimuth)
    grad_B = np.tan(np.radians(dip))
    Z = -grad_B * (np.sin(azir)*X + np.cos(azir)*Y)
    #
    # Affine Transformation
    transform = Affine.translation(x_UpperLeft, y_UpperLeft) * Affine.scale(dxy, -dxy) # affine transformation matrix
    #
    bbox = 'left={:.0f}__bottom={:.0f}__right={:.0f}__top={:.0f}'.format(x.min(), y.min(), x.max(), y.max())
    file = 'inclined_plane_'+'dip={}_azi={}_'.format(dip, azimuth)+bbox+'.tif'
    fullpath = path_to_dir + file
    #
    new_dataset = rasterio.open(
        fullpath,
        'w',
        driver='GTiff',
        height=Z.shape[0],
        width=Z.shape[1],
        count=1,
        dtype=Z.dtype,
        crs='+proj=latlong',
        transform=transform,
        )
    #
    new_dataset.write(Z, 1)
    new_dataset.close()


def make_cone(dip, nx, ny, dxy, path_to_dir):
    # Make cone dem
    # dip: how steep (degrees down from horizontal)
    # lat,lon will be equal to x,y
    #
    x_UpperLeft = -(nx//2)*dxy
    y_UpperLeft = (ny//2)*dxy
    x = x_UpperLeft + np.arange(nx)*dxy
    y = y_UpperLeft - np.arange(ny)*dxy
    X,Y = np.meshgrid(x,y)
    R = (X**2 + Y**2)**0.5
    #
    grad_B = np.tan(np.radians(dip))
    Z = -grad_B * R
    #
    # Affine Transformation
    transform = Affine.translation(x_UpperLeft, y_UpperLeft) * Affine.scale(dxy, -dxy) # affine transformation matrix
    #
    bbox = 'left={:.0f}__bottom={:.0f}__right={:.0f}__top={:.0f}'.format(x.min(), y.min(), x.max(), y.max())
    file = 'cone_'+'dip={}_'.format(dip)+bbox+'.tif'
    fullpath = path_to_dir + file
    #
    new_dataset = rasterio.open(
        fullpath,
        'w',
        driver='GTiff',
        height=Z.shape[0],
        width=Z.shape[1],
        count=1,
        dtype=Z.dtype,
        crs='+proj=latlong',
        transform=transform,
        )
    #
    new_dataset.write(Z, 1)
    new_dataset.close()


def smooth_topo(n_times):
    for i in range(n_times):
        print('| --- Smoothing Topography --- |')
        g.B0 = m.local_mean(g.B0)



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
    fully_crust = phi_true_crust >= 1.0 # if all neighborhood cells true
    core_locked = cryst_core >= p.phi_max # if all neighborhood cells true
    stalled = m.local_and_EWNS(jeffreys_efficiency == 0) # EWNS neighbors all stopped
    freezeout = np.logical_and(stalled, np.logical_or(core_locked, fully_crust))
    #
    h_to_freeze = p.fraction_to_freeze * h[freezeout]
    B_n[freezeout] += h_to_freeze
    h[freezeout] -= h_to_freeze
    dBdt[freezeout] = 0
    #
    #ti[freezeout] = t_n
    # R, Rs, abs_U_s all zero already in freezout region
    # no return needed







#

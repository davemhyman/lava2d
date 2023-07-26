import numpy as np
import scipy.interpolate as sci
import pandas as pd
import glob
import os

from globals import params as p

################################################################################
#
################################################################################


def spatial_density(x,y,x0,y0,x1,y1,W):
    Dx = x1-x0
    Dy = y1-y0
    L = np.sqrt(Dx**2 + Dy**2)
    # Adjust Dx for theta
    if Dx == 0:
        Dx = p.pos_eps
    theta = np.arctan(Dy/Dx)
    cos_th = np.cos(theta)
    sin_th = np.sin(theta)
    xn = x - 0.5*(x1+x0)
    yn = y - 0.5*(y1+y0)
    u = xn*cos_th + yn*sin_th
    v = -xn*sin_th + yn*cos_th
    #
    # Adjust L,W so f>0 (only if f has compact support)
    L = max(L, 1.01*min(p.dx,p.dy))
    W = max(W, 1.01*min(p.dx,p.dy))
    #
    f = (np.sign(u+0.5*L) - np.sign(u-0.5*L)) * np.maximum(1-(v/(0.5*W))**2, 0)
    f_tot = np.sum(f)
    if f_tot == 0:
        return f # since its already zero
    else:
        return f / (np.sum(f)*p.dx*p.dy) # required for correct total (f.sum()*dx*dy = 1)



def read_source_data():
    p.vent_param_splines=[]
    files_to_search = os.path.join(p.path_to_vent_files, 'vent_*.txt')
    files = glob.glob(files_to_search)
    for n in range(len(files)):
        # open each vent file:
        db = pd.read_csv(files[n], header = 0, delimiter = '\t')
        db['discharge'] = db['discharge'] / (1.0-p.porosity) # convert from DRE to Bulk Effusion Rate
        db.dropna(inplace = True)
        data = db.to_numpy().T
        times = data[0]
        params = data[1:]
        # construct linear interpolant for each fissure parameter
        Pn = sci.interp1d(times, params, kind = 'linear', bounds_error = False, fill_value = (params[:,0], params[:,-1]))
        #
        p.vent_param_splines.append(Pn)


def read_list_data(vents_list):
    p.vent_param_splines = []
    for vent in vents_list:
        db = pd.DataFrame(vent, columns=['time', 'x0', 'y0', 'x1', 'y1', 'width', 'discharge'])
        db['discharge'] = db['discharge'] / (1.0-p.porosity) # convert from DRE to Bulk Effusion Rate
        db.dropna(inplace=True)
        data = db.to_numpy().T
        times = data[0]
        params = data[1:]
        # construct linear interpolant for each fissure parameter
        Pn = sci.interp1d(times, params, kind='linear', bounds_error=False, fill_value=(params[:,0], params[:,-1]))
        #
        p.vent_param_splines.append(Pn)



def source_term(x,y,t):
    src = np.zeros(x.shape)
    I = np.full(x.shape, False) # boolean array for vents
    x_UpperLeft, y_UpperLeft = x[0,0], y[0,0]
    ny, nx = x.shape
    #
    I = src > 0 # boolean array for vents
    dxy_min = min(p.dx,p.dy)
    #
    n_vents = len(p.vent_param_splines)
    for n in range(n_vents):
        # at some time t:
        x0_n, y0_n, x1_n, y1_n, W_n, Q_n = p.vent_param_splines[n](t) # scalars
        if Q_n != 0:
            L_n = np.sqrt((x1_n-x0_n)**2 + (y1_n-y0_n)**2)
            bw_n = int(max(L_n,W_n) / dxy_min) + 1
            xavg_n = 0.5*(x0_n+x1_n)
            yavg_n = 0.5*(y0_n+y1_n)
            id_xavg_n = int((xavg_n - x_UpperLeft) / p.dx)
            id_yavg_n = int((y_UpperLeft - yavg_n) / p.dy)
            #
            r_top = max(id_yavg_n - bw_n, 0)
            r_bottom = min(id_yavg_n + bw_n+1, ny)
            c_left = max(id_xavg_n - bw_n, 0)
            c_right = min(id_xavg_n + bw_n+1, nx)
            #
            I[r_top:r_bottom, c_left:c_right] = True
            #
            # only compute source term over vent region n
            src[I] += Q_n * spatial_density(x[I], y[I], x0_n, y0_n, x1_n, y1_n, W_n)
        #
    return src











#

import numpy as np
import scipy.interpolate as sci
import pandas as pd
import glob
import numexpr as ne

from globals import params as p
from globals import grids as g
import thermal as therm
import topo
import rheo

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
    #u = ne.evaluate('xn*cos_th + yn*sin_th')
    #v = ne.evaluate('-xn*sin_th + yn*cos_th')
    u = xn*cos_th + yn*sin_th
    v = -xn*sin_th + yn*cos_th
    #
    # Adjust L,W so f>0 (only if f has compact support)
    L = max(L, 1.01*min(p.dx,p.dy))
    W = max(W, 1.01*min(p.dx,p.dy))
    #
    #f = 0.75/(L*W) * (np.sign(u+0.5*L) - np.sign(u-0.5*L)) * np.maximum(1-(v/(0.5*W))**2,0)
    f = (np.sign(u+0.5*L) - np.sign(u-0.5*L)) * np.maximum(1-(v/(0.5*W))**2, 0)
    return f / (np.sum(f)*p.dx*p.dy + p.pos_eps) # required for correct total (f.sum()*dx*dy = 1)



def read_source_data():
    p.vent_param_splines=[]
    files = glob.glob(p.path_to_vent_files + 'vent_*.txt')
    for n in range(len(files)):
        # open each vent file:
        db = pd.read_csv(files[n], header = 0, delimiter = '\t')
        db.dropna(inplace = True)
        data = db.to_numpy().T
        times = data[0]
        params = data[1:]
        # construct linear interpolant for each fissure parameter
        Pn = sci.interp1d(times, params, kind = 'linear', bounds_error = False, fill_value = (params[:,0], params[:,-1]))
        #
        p.vent_param_splines.append(Pn)





def source_term(x,y,t,I):
    # at some time t:
    x_active,y_active = x[I], y[I] # (n_active,)
    #
    n_vents = len(p.vent_param_splines)
    q = np.zeros(x_active.shape) # (n_active,)
    Q_tot = 0
    for n in range(n_vents):
        x0,y0,x1,y1,W,Qn = p.vent_param_splines[n](t) # scalars
        fn = spatial_density(x_active,y_active,x0,y0,x1,y1,W) # (n_active,)
        q += Qn * fn # (n_active,)
        Q_tot += Qn
    # Need to correct to ensure np.sum(q)*p.dx*p.dy = Q_tot
    # as a safety from overlaps
    src = np.zeros(x.shape)
    src[I] = Q_tot * q / (np.sum(q)*p.dx*p.dy + p.pos_eps) # (rows, cols)
    return src


















#

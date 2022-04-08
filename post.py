import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import colorbar as cbar
import matplotlib.colors as c
#from mg_fas import *
import rheo
import thermal as therm
import vents

from globals import params as p
from globals import grids as g


# Post Processing


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = c.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n))
        )
    return new_cmap






def make_flow_maps(fields, dz = 10):
    #
    B_contours = np.arange(np.floor(g.B0.min()/dz), np.ceil(g.B0.max()/dz) + 2) * dz
    new_cmap = truncate_colormap(plt.get_cmap('terrain'), 0.25, 0.9)
    I = g.h_n > p.tiny_flow
    #
    x_vent = g.x[vents.source_term(g.x,g.y,p.t_max,np.isfinite(g.x))>p.pos_eps]
    y_vent = g.y[vents.source_term(g.x,g.y,p.t_max,np.isfinite(g.x))>p.pos_eps]
    #
    for type in fields:
        if type == 'stoch_bed':
            #
            fig, ax1 = plt.subplots(figsize=(16,8))
            plt.contourf(g.x, g.y, g.B_n-g.B0, 51, cmap= 'terrain')
            plt.axis('equal')
            plt.plot(x_vent, y_vent, 'r^')
            plt.title('Stochastic Bed Variation (m)')
            plt.colorbar();plt.show()
            #
        elif type == 'height':
            #
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            cB = plt.contour(g.x, g.y, g.B_n, B_contours, linewidths = .2, colors = 'k', zorder = 2)
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.h_n, ~I), cmap = 'magma', shading = 'auto', zorder = 3);plt.colorbar()
            cid = plt.contour(g.x, g.y, g.h_n, 21, linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('Flow Thickness (m)')
            plt.show()
            #
        elif type == 'surface_height':
            #
            S_I = I*(g.B_n+g.h_n)
            #
            S_contours = np.arange(np.floor(np.min(S_I[I])), np.ceil(np.max(S_I)), 1)
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.h_n, ~I), cmap = 'magma', shading = 'auto', zorder = 2);plt.colorbar()
            cid = plt.contour(g.x, g.y, S_I, S_contours, linewidths = .1, colors = 'k', zorder = 3)
            cB = plt.contour(g.x, g.y, g.B_n+g.h_n, B_contours, linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('Flow Thickness (m)')
            plt.show()
            #
        elif type == 'log-height':
            #
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            cB = plt.contour(g.x, g.y, g.B_n, B_contours, linewidths = .2, colors = 'k', zorder = 2)
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.h_n,~I), cmap = 'magma', shading = 'auto', norm = c.LogNorm(), vmin = p.tiny_flow, vmax = np.ceil(g.h_n.max()), zorder = 3);plt.colorbar()
            cid = plt.contour(g.x, g.y, g.h_n, p.pos_eps+np.arange(0, np.ceil(g.h_n.max()), 0.1), linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('Flow Thickness (m)')
            plt.show()
            #
        elif type == 't_inundation':
            #
            tmax = g.t_inundation.max()
            if tmax/86400 < 2:
                ti_contours = np.arange(0,tmax/3600 + 1, 1) # 1 hr intervals
            elif tmax/86400 >= 2 and tmax/86400 < 4:
                ti_contours = np.arange(0,tmax/3600 + 2, 2) # 2 hr intervals
            elif tmax/86400 >= 4 and tmax/86400 < 8:
                ti_contours = np.arange(0,tmax/3600 + 6, 6) # 6 hr intervals
            elif tmax/86400 >= 8 and tmax/86400 < 15:
                ti_contours = np.arange(0,tmax/3600 + 12, 12) # 12 hr intervals
            elif tmax/86400 >= 30:
                ti_contours = np.arange(0,tmax/3600 + 24, 24) # 24 hr intervals
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            cB = plt.contour(g.x, g.y, g.B_n, B_contours, linewidths = .2, colors = 'k', zorder = 2)
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.t_inundation/3600,~I), cmap = 'magma', shading = 'auto', zorder = 3);plt.colorbar()
            cid = plt.contour(g.x, g.y, g.t_inundation/3600, ti_contours, linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('Inundation Time (hrs)')
            plt.show()
            #
        elif type == 't_erupted':
            #
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            cB = plt.contour(g.x, g.y, g.B_n, B_contours, linewidths = .2, colors = 'k', zorder = 2)
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.t_erupted/60,~I), cmap = 'magma_r', shading = 'auto', zorder = 3);plt.colorbar()
            cid = plt.contour(g.x, g.y, g.t_erupted/60, np.arange(0,np.ceil(g.t_erupted.max()/60), 5), linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('Inundation Time (min)')
            plt.show()
            #
        elif type == 'surf_T':
            #
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            cB = plt.contour(g.x, g.y, g.B_n, B_contours, linewidths = .2, colors = 'k', zorder = 2)
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.surface_T_n,~I), vmin = 500, vmax = 1400, cmap = 'hot', shading = 'auto', zorder = 3);plt.colorbar()
            cid = plt.contour(g.x, g.y, g.surface_T_n, np.arange(500, 1500,10), linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('Surface Temperature (K)')
            plt.show()
            #
        elif type == 'viscosity':
            #
            vmin = 10**(np.floor(np.log10(g.mu_core_n.min())))
            vmax = 1e4 * vmin
            log_contours = (10.**np.arange(np.log10(vmin),np.log10(vmax))[:,None] * np.arange(1,10)[None,:]).reshape(-1)
            log_contours = np.hstack([log_contours, vmax])
            #
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            cB = plt.contour(g.x, g.y, g.B_n, B_contours, linewidths = .2, colors = 'k', zorder = 2)
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.mu_core_n,~I), norm = c.LogNorm(), vmin =vmin, vmax = vmax, cmap = 'hot_r', shading = 'auto', zorder = 3);plt.colorbar()
            cid = plt.contour(g.x, g.y, g.mu_core_n, log_contours, linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('viscosity (Pa s)')
            plt.show()
            #
        elif type == 'Re':
            #
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            cB = plt.contour(g.x, g.y, g.B_n, B_contours, linewidths = .2, colors = 'k', zorder = 2)
            vmin = 1e-3
            vmax = 10**(np.ceil(np.log10(g.Re_n.max())))
            log_contours = (10.**np.arange(np.log10(vmin),np.log10(vmax))[:,None] * np.arange(1,10)[None,:]).reshape(-1)
            log_contours = np.hstack([log_contours, vmax])
            #
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.Re_n,~I), cmap = 'viridis', shading = 'auto', norm = c.LogNorm(), vmin = vmin, vmax = vmax, zorder = 3);plt.colorbar()
            cid = plt.contour(g.x, g.y, g.Re_n, log_contours, linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('Reynolds Number')
            plt.show()
        elif type == 'Bn':
            #
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            cB = plt.contour(g.x, g.y, g.B_n, B_contours, linewidths = .2, colors = 'k', zorder = 2)
            vmin = 1e-4
            vmax = 1
            log_contours = (10.**np.arange(np.log10(vmin),np.log10(vmax))[:,None] * np.arange(1,10)[None,:]).reshape(-1)
            log_contours = np.hstack([log_contours, vmax])
            #
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.Bn_core_n,~I), cmap = 'viridis', shading = 'auto', norm = c.LogNorm(), vmin = vmin, vmax = vmax, zorder = 3);plt.colorbar()
            cid = plt.contour(g.x, g.y, g.Bn_core_n, log_contours, linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('Bingham Number')
            plt.show()
            #
        elif type == 'K':
            #
            fig, ax1 = plt.subplots(figsize=(16,8))
            cfB = plt.contourf(g.x, g.y, g.B_n, B_contours, cmap= new_cmap, zorder = 1)
            cB = plt.contour(g.x, g.y, g.B_n, B_contours, linewidths = .2, colors = 'k', zorder = 2)
            vmin = 1e-2
            vmax = 10**(np.ceil(np.log10(g.K_n.max())))
            log_contours = (10.**np.arange(np.log10(vmin),np.log10(vmax))[:,None] * np.arange(1,10)[None,:]).reshape(-1)
            log_contours = np.hstack([log_contours, vmax])
            #
            pid = plt.pcolormesh(g.x, g.y, np.ma.masked_array(g.K_n,~I), cmap = 'viridis', shading = 'auto', norm = c.LogNorm(), vmin = vmin, vmax = vmax, zorder = 3);plt.colorbar()
            cid = plt.contour(g.x, g.y, g.K_n, log_contours, linewidths = .2, colors = 'k', zorder = 4)
            axeq = plt.axis('equal')
            U = np.sqrt(g.ux_n**2 + g.uy_n**2)
            qvr = ax1.quiver(g.x[::2,::2],g.y[::2,::2],g.ux_n[::2,::2], g.uy_n[::2,::2], scale_units='x', scale = U.max()/(10*p.dx), units = 'x', width = 0.1*p.dx, color='blue', zorder = 5)
            vid = plt.plot(x_vent, y_vent, 'r^', zorder = 6)
            tid = plt.title('Diffusivity (m2 s-1)')
            plt.show()










def bingham_number(tau_0):
    S = g.B_n + g.h_n
    beta = tau_0 / (p.lava_density * p.g * g.h_n + p.pos_eps)
    beta_ij = beta[1:-1,1:-1] # center
    S_ij = S[1:-1,1:-1] # center
    #
    dSdx_R = (S[1:-1,2:] - S_ij) / p.dx # x-deriv on right
    dSdx_L = (S_ij - S[1:-1,0:-2]) / p.dx # x-deriv on left
    dSdy_U = (S[0:-2,1:-1] - S_ij) / p.dy # y-deriv on upper
    dSdy_D = (S_ij - S[2:,1:-1]) / p.dy # y-deriv on lower (down)
    dSdx_UD = 0.5 * (dSdx_R + dSdx_L) # x-deriv on upper/lower
    dSdy_LR = 0.5 * (dSdy_U + dSdy_D) # x-deriv on upper/lower
    grad_S_R = np.sqrt(dSdx_R**2 + dSdy_LR**2) + p.pos_eps # div-by-zero protection
    grad_S_L = np.sqrt(dSdx_L**2 + dSdy_LR**2) + p.pos_eps
    grad_S_U = np.sqrt(dSdx_UD**2 + dSdy_U**2) + p.pos_eps
    grad_S_D = np.sqrt(dSdx_UD**2 + dSdy_D**2) + p.pos_eps
    #
    Bn_R = 0.5 * (beta_ij/grad_S_R + beta[1:-1,2:]/grad_S_R)
    Bn_L = 0.5 * (beta_ij/grad_S_L + beta[1:-1,0:-2]/grad_S_L)
    Bn_U = 0.5 * (beta_ij/grad_S_U + beta[0:-2,1:-1]/grad_S_U)
    Bn_D = 0.5 * (beta_ij/grad_S_D + beta[2:,1:-1]/grad_S_D)
    #
    out = np.zeros(g.h_n.shape)
    out[1:-1,1:-1] = 0.25 * (Bn_R + Bn_L + Bn_U + Bn_D)
    return out






def velocity2d(h_n, B_n, R_n):
    #
    S = B_n + h_n
    Q = R_n * h_n**2
    #
    h_ij = h_n[1:-1,1:-1] # center
    S_ij = S[1:-1,1:-1] # center
    Q_ij = Q[1:-1,1:-1] # center
    #
    dSdx_R = (S[1:-1,2:] - S_ij) / p.dx # x-deriv on right
    dSdx_L = (S_ij - S[1:-1,0:-2]) / p.dx # x-deriv on left
    dSdy_U = (S[0:-2,1:-1] - S_ij) / p.dy # y-deriv on upper
    dSdy_D = (S_ij - S[2:,1:-1]) / p.dy # y-deriv on lower (down)
    #
    ux = np.zeros(h_n.shape)
    uy = np.zeros(h_n.shape)
    ux[1:-1,1:-1] = -((Q_ij + Q[1:-1,2:]) * dSdx_R + (Q_ij + Q[1:-1,0:-2]) * dSdx_L) / 4.
    uy[1:-1,1:-1] = -((Q_ij + Q[0:-2,1:-1]) * dSdy_U + (Q_ij + Q[2:,1:-1]) * dSdy_D) / 4.
    #
    return ux, uy


def surface_velocity2d(h_n, B_n, Rs_n):
    #us_x,us_y = surface_velocity2d(h_next, B_next, Rs_next)
    #
    S = B_n + h_n
    Q = Rs_n*h_n**2
    #
    h_ij = h_n[1:-1,1:-1] # center
    S_ij = S[1:-1,1:-1] # center
    Q_ij = Q[1:-1,1:-1] # center
    #
    dSdx_E = (S[1:-1,2:] - S_ij) / p.dx # x-deriv on right
    dSdx_W = (S_ij - S[1:-1,0:-2]) / p.dx # x-deriv on left
    dSdy_N = (S[0:-2,1:-1] - S_ij) / p.dy # y-deriv on upper
    dSdy_S = (S_ij - S[2:,1:-1]) / p.dy # y-deriv on lower (down)
    #
    ux = np.zeros(h_n.shape)
    uy = np.zeros(h_n.shape)
    ux[1:-1,1:-1] = -((Q_ij + Q[1:-1,2:]) * dSdx_E + (Q_ij + Q[1:-1,0:-2]) * dSdx_W) / 4.
    uy[1:-1,1:-1] = -((Q_ij + Q[0:-2,1:-1]) * dSdy_N + (Q_ij + Q[2:,1:-1]) * dSdy_S) / 4.
    #
    return ux, uy


def make_all_physics_grids(t_n):
    #
    I = np.isfinite(g.x)
    cryst_core = rheo.cryst_avrami(t_n, g.t_erupted, I)
    g.mu_core_n = rheo.viscosity(p.core_temperature, cryst_core)
    g.tau_0_core_n = rheo.yield_stress(cryst_core)
    g.Bn_core_n = bingham_number(g.tau_0_core_n)
    g.phi_S = therm.surface_BL(t_n, g.t_erupted, g.h_n, p.core_temperature)
    g.R_n, g.Rs_n, fluidity_n, abs_grad_S_n = rheo.rheo_factor_bl_S(g.h_n, g.B_n, g.phi_S, cryst_core)
    kinematic_viscosity =  np.sqrt(2) * g.mu_core_n / p.lava_density
    g.Pr_n = kinematic_viscosity / p.lava_diffusivity
    g.Pe_n = g.abs_Usurf * g.h_n / p.lava_diffusivity
    g.Fo_n = p.lava_diffusivity * (t_n - g.t_inundation) / np.maximum(g.h_n, p.pos_eps)**2 + p.pos_eps
    #
    g.K_n = g.R_n * g.h_n**3
    g.ux_n, g.uy_n = velocity2d(g.h_n, g.B_n, g.R_n)
    g.Re_n = 3 * g.R_n * np.sqrt(g.ux_n**2 + g.uy_n**2) * g.h_n / p.g
    #
    valid = (g.h_n > p.tiny_flow)
    tv_s = t_n - g.t_erupted[valid] + p.pos_eps # (n_valid_x)
    g.surface_T_n = p.atm_temperature + np.zeros(g.h_n.shape, dtype = g.h_n.dtype) # initialize output space
    if np.any(valid):
        Ts, hc, T_inf, n_iter = therm.solve_T_surf_robin_problem(tv_s, p.core_temperature)
        g.surface_T_n[valid] = Ts


#

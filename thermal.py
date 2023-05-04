import numpy as np
from scipy.special import erfc

from globals import params as p
from globals import grids as g
import rheo
################################################################################
# Energy Equation Components
################################################################################
def surf_heat_flux(T):
    sb_eff = p.lava_emissivity * p.stefan_boltzmann
    d2_rad = 12 * sb_eff * T**2
    d_rad = d2_rad * T / 3.
    rad = d_rad * T / 4 - sb_eff * p.atm_temperature**4
    #
    conv = p.h_conv * (T - p.atm_temperature)
    #
    q = rad + conv
    dq = d_rad + p.h_conv
    #
    return q, dq, d2_rad

def erfinv(x):
    a = 0.147
    b = 2./(np.pi * a)
    ln = np.log(1.0 - x**2)
    q1 = b + 0.5 * ln
    q2 = ln / a
    sqrt1 = np.sqrt(q1**2 - q2)
    sqrt2 = np.sqrt(sqrt1 - q1)
    return np.sign(x) * sqrt2

def Theta(eta, xi):
    tau = 1. / (1+0.5*(eta + xi))
    P_tau = -1.26551223 + tau*(1.00002368 + tau*(.37409196 + tau*(.09678418 + tau*(-.18628806 + tau*(.27886807 + tau*(-1.13520398 + tau*(1.48851587 + tau*(-.82215223 + tau*.17087277))))))))
    return erfc(eta) - tau*np.exp(-eta**2 + P_tau)

def solve_T_surf_robin_problem(t_s, T_c, tol = 1e-6):
    #
    Lk = np.sqrt(t_s / (p.lava_density * p.lava_specific_heat * p.lava_conductivity))
    # initial guess
    sb_eff = p.lava_emissivity * p.stefan_boltzmann
    dq_core = 4 * sb_eff * T_c**3 + p.h_conv
    dq_atm = 4 * sb_eff * p.atm_temperature**3 + p.h_conv
    Ts = p.atm_temperature - (T_c - p.atm_temperature) / (Theta(0,dq_core*Lk) / (Theta(0,dq_atm*Lk)-1) - 1)
    #
    q, dq, d2q = surf_heat_flux(Ts) # surface heat flux and derivatives
    xi = dq * Lk
    Theta_surf = Theta(0,xi)
    f = Theta_surf - (Ts - T_c) / (Ts - T_c - q/dq)
    #
    abs_err = abs(f).max()
    max_it = 10
    n = 0
    while (abs_err >= tol) and (n <= max_it):
        # Newton's Method
        # given current vals: Ts, q, dq, xi, Theta_surf, f
        # Build dTs
        df = 2*d2q*Lk/np.sqrt(np.pi) * (1 - np.sqrt(np.pi)*xi*(1-Theta_surf)) - ((Ts-T_c)*(dq**2 - q*d2q) - q*dq) / ((Ts-T_c)*dq - q)**2
        dTs = - f / df
        # update
        Ts = np.clip(Ts + dTs, p.atm_temperature, T_c)
        #
        # compute intermediate vals for next iter:
        q, dq, d2q = surf_heat_flux(Ts)
        xi = dq * Lk
        Theta_surf = Theta(0,xi)
        f = Theta_surf - (Ts - T_c) / (Ts - T_c - q/dq)
        #
        abs_err = abs(f).max()
        n+=1
    #
    if n >= max_it:
        print('Failed to Converge by max_it')
    #
    # returns:
    # surface temp, local effective h_conv, local effective T_inf, n_iterations
    return Ts, dq, Ts - q / dq, n


def surface_BL_simple(t, te, h, T_core):
    # t_s: time since eruption of surface material for each cell
    # h: gravity current height at each cell
    phi_S = np.zeros(h.shape, dtype = h.dtype) # initialize output space # # all h < p.tiny_flow fully solid
    #
    valid = (h > p.tiny_flow)
    #
    if np.any(valid):
        #
        tv_s = t - te[valid] # (n_valid)
        hv = h[valid] # (n_valid)
        #
        # SURFACE BL
        #
        T_stiff = rheo.T_stiffened(T_core) # temp at which visc = 2*visc_core
        eta_inf = erfinv((T_stiff-p.atm_temperature)/(T_core-p.atm_temperature))
        phi_S[valid] = eta_inf * np.sqrt(4 * p.lava_diffusivity * tv_s) / hv
        #
    #
    return phi_S


#

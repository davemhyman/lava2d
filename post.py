import numpy as np

from globals import params as p
from globals import grids as g
import rheo
import thermal as therm

# Post Processing

def make_all_physics_grids(t_n):
    #
    I = np.isfinite(g.x)
    cryst_core = rheo.cryst_avrami(t_n, g.t_erupted, I)
    g.mu_core_n = rheo.viscosity(p.core_temperature, cryst_core)
    g.tau_0_core_n = rheo.yield_stress(cryst_core)
    g.Bn_core_n = bingham_number(g.tau_0_core_n)
    g.phi_S = therm.surface_BL_simple(t_n, g.t_erupted, g.h_n, p.core_temperature)
    g.R_n, g.Rs_n, fluidity_n, abs_grad_S_n, jeffreys_efficiency = rheo.rheo_factor_bl_S_all_outputs(g.h_n, g.B_n, g.phi_S, cryst_core)
    kinematic_viscosity =  np.sqrt(2) * g.mu_core_n / p.lava_density
    g.Pr_n = kinematic_viscosity / p.lava_diffusivity
    g.Pe_n = g.abs_Usurf * g.h_n / p.lava_diffusivity
    g.Fo_n = p.lava_diffusivity * (t_n - g.t_inundation) / np.maximum(g.h_n, p.pos_eps)**2 + p.pos_eps
    #
    g.K_n = g.R_n * g.h_n**3
    g.ux_n, g.uy_n = velocity2d(g.h_n, g.B_n, g.R_n)
    g.Re_n = 3 * g.R_n * np.sqrt(g.ux_n**2 + g.uy_n**2) * g.h_n / p.g
    #
    # Compute surface T for unfrozen and frozen
    dB = g.B_n - g.B0
    valid = (dB + g.h_n > p.tiny_flow)
    tv_s = t_n - g.t_erupted[valid] + p.pos_eps # (n_valid_x)
    g.surface_T_n = p.atm_temperature + np.zeros(g.h_n.shape, dtype = g.h_n.dtype) # initialize output space
    if np.any(valid):
        Ts, hc, T_inf, n_iter = therm.solve_T_surf_robin_problem(tv_s, p.core_temperature)
        g.surface_T_n[valid] = Ts


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

#

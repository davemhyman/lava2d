import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import interpolate
import numexpr as ne

from globals import params as p
from globals import grids as g
import thermal as therm
import vents
import topo
import rheo

################################################################################
# Time Updating Operations
################################################################################
def update_inundation_times(t_inundation, h_n, t, dt):
    # most recently inundated cells should be assigned an inundation time of
    #     t_inundation = t - 0.5*dt (t-t_inundation = 0.5*dt passed to T_avg).
    # All un-inundated cells should be assigned an inundation time of
    #     t_inundation = t (t-t_inundation = 0 passed to thermal problems)
    #
    not_inundated_last = np.isclose(t_inundation, t-dt, rtol = 0.0, atol = .001*dt) # un-inundated cells last time step
    inundated_now = h_n > p.tiny_flow # inundated cells
    update_t_in =  np.logical_and(not_inundated_last, inundated_now) # inundated during this time step
    #
    t_inundation[not_inundated_last] = t # all un-inundated cells --> time t
    t_inundation[update_t_in] = t - 0.5*dt # all recently inundated cells --> time t-dt/2
    return t_inundation





def update_te_from_advection(te, Rsurf_n, neg_grad_S, h2_RRUU, h2_LLDD, max_h2, q_n, dxy, t_n, dt, is_uphill_flow, I, I_ij):
    # Implicitly global:
    # Rsurf_n, neg_grad_S, h2_RRUU, h2_LLDD, max_h2, q_n, dxy, t_n, dt, is_uphill_flow, I, I_ij
    #
    # global te
    #
    # t_erupted
    te_ij = te[I] # (n_active,) # Q[I] == Q[2:-2,2:-2][I_ij]
    Rsurf_ij = Rsurf_n[I] # (n_active,)
    Rsurf_RLUD = np.array([Rsurf_n[2:-2,3:-1][I_ij], Rsurf_n[2:-2,1:-3][I_ij], Rsurf_n[1:-3,2:-2][I_ij], Rsurf_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    Rsurf = 0.5 * (Rsurf_ij[None,:] + Rsurf_RLUD) # (4, n_active) # order: E, W, N, S
    # t_erupted MUSCL Reconstruction:
    all_te_states = np.array(muscl_linear(te, I, flux_lim_superbee)) # (8, n_active)
    te_RRUU = all_te_states[0::2] # (4, n_active)
    te_LLDD = all_te_states[1::2] # (4, n_active)
    # --------------------------------------------------------------------------
    # Surface Velocities:
    neg_Rsurf_grad_S = Rsurf * neg_grad_S
    Usurf = neg_Rsurf_grad_S * 0.5 * (h2_RRUU + h2_LLDD) # (4, n_active) # order: E, W, N, S
    asurf = abs(neg_Rsurf_grad_S) * max_h2 # (4, n_active) # order: E, W, N, S
    abs_Usurf = np.zeros(te.shape, dtype = te.dtype)
    abs_Usurf[I] = np.max(abs(Usurf), axis = 0)
    # --------------------------------------------------------------------------
    # Update t_erupted
    all_fluxes_te = 0.5 * (
        (neg_Rsurf_grad_S*h2_RRUU - asurf)*te_RRUU +
        (neg_Rsurf_grad_S*h2_LLDD + asurf)*te_LLDD
        ) # (4, n_active)
    #
    all_fluxes_te[is_uphill_flow] = 0 # (4, n_active) # order: E, W, N, S
    Usurf[is_uphill_flow] = 0 # (4, n_active) # order: E, W, N, S
    #
    #neg_advection_te = te_ij*np.sum(Usurf/dxy[:,None], axis = 0) - np.sum(all_fluxes/dxy[:,None], axis = 0)
    neg_advection_te = np.sum((te_ij[None,:]*Usurf - all_fluxes_te) / dxy[:,None], axis = 0) # numerically the same as commented version
    te[I] += dt * neg_advection_te
    te[q_n > 0] = t_n + dt # ensure source conditon  # (rows, cols)
    #
    return abs_Usurf



def update_te_from_inundation(te, Rsurf_n, abs_grad_S_n, h, R_n, ti, t_n, dt, I):
    # Implicitly global:
    # Rsurf_n, abs_grad_S_n, h, R_n, ti, t_n, dt, I
    #
    # global te
    #
    # --------------------------------------------------------------------------
    # Surface Speed
    abs_U_s = np.zeros(h.shape, dtype = h.dtype)
    abs_U_s[I] = Rsurf_n[I] * abs_grad_S_n[I] * h[I]**2
    #
    # --------------------------------------------------------------------------
    # Update t_erupted
    speed_ratio_ij = local_mean((R_n + p.pos_eps) / (Rsurf_n + p.pos_eps))[I]
    te_ij_update = t_n - speed_ratio_ij * ti[I] # assume surface travel time is ~ 2/3 of bulk current
    te_ij_update[q_ij > 0] = t_n + dt # ensure source conditon  # (rows, cols)
    te[I] = te_ij_update
    #
    return abs_Usurf



################################################################################
# Adaptive Grid Subset Tools
################################################################################
def subset_bounds(t, buffer_cells = 5):
    #
    ny_global, nx_global = g.h_n.shape
    #
    x_UpperLeft, y_UpperLeft = g.x[0,0], g.y[0,0]
    #
    I = g.h_n > 0
    #
    '''
    # Correct for vent areas
    I[180:225,65:95] = True
    #
    # need faster vent area corrector
    # static corrector is pretty fast if vents are close together
    # cannot compute corrector region every time -- too slow!
    # could make a global grid g.t_on with vent area cells having the value of time when they turn on
    # similar grid for when they shut off (g.t_off)
    # (g.t_on <= t_n) and (t_n <= g.t_off) demarcates vent areas where discharge is ongoing
    # g.discharge_ongoing = np.logical_and(g.t_on <= t_n, g.t_off >= t_n)
    # I = np.logical_or(g.h_n > 0, g.discharge_ongoing)
    #
    #
    #
    #
    n_vents = len(p.vent_param_splines)
    dxy_min = min(p.dx,p.dy)
    for n in range(n_vents):
        x0,y0,x1,y1,W,Q_n = p.vent_param_splines[n](t) # scalars
        if Q_n != 0:
            L = np.sqrt((x1-x0)**2 + (y1-y0)**2)
            bw = int(0.5 * max(L,W) / dxy_min) + 1
            xavg = 0.5*(x0+x1)
            yavg = 0.5*(y0+y1)
            #
            id_xavg = int((xavg - x_UpperLeft) / p.dx)
            id_yavg = int((y_UpperLeft - yavg) / p.dy)
            #
            r_top = max(id_yavg - bw, 0)
            r_bottom = min(id_yavg + bw+1, ny_global)
            c_left = max(id_xavg - bw, 0)
            c_right = min(id_xavg + bw+1, nx_global)
            #
            I[r_top:r_bottom, c_left:c_right] = True
        #
    #
    '''
    #
    if ~np.any(I):
        # h = 0, q = 0 ==> nothing will happen
        # --> make active at 1 cell in the grid center
        sh = g.x.shape
        I[sh[0]//2, sh[1]//2] = True
    #
    active = I.astype(g.h_n.dtype)
    #
    idx_x = np.where(np.sum(active, axis = 0) > 0)[0]
    idx_y = np.where(np.sum(active, axis = -1) > 0)[0]
    #
    col_min = idx_x.min() - buffer_cells
    row_min = idx_y.min() - buffer_cells
    col_max = idx_x.max() + buffer_cells
    row_max = idx_y.max() + buffer_cells
    nx_a = col_max - col_min
    ny_a = row_max - row_min
    nx = int(2**np.ceil(np.log2(nx_a))) + 1 # odd
    ny = int(2**np.ceil(np.log2(ny_a))) + 1 # odd
    wx = (nx-1)//2 # even
    wy = (ny-1)//2 # even
    col_0 = (col_min + col_max)//2
    row_0 = (row_min + row_max)//2
    #
    r_start = row_0 - wy
    r_stop = row_0 + wy + 1
    c_start = col_0 - wx
    c_stop = col_0 + wx + 1
    #
    # Corrections if applicable
    if r_start < 0:
        r_start = 0
        r_stop = min(ny, ny_global) # both are 2**n + 1
    if c_start < 0:
        c_start = 0
        c_stop = min(nx, nx_global) # both are 2**n + 1
    if r_stop > ny_global:
        r_stop = ny_global
        r_start = max(ny_global-ny, 0)
    if c_stop > nx_global:
        c_stop = nx_global
        c_start = max(nx_global-nx, 0)
    #
    return r_start, r_stop, c_start, c_stop

def subset_grids(bounds):
    #
    r_start, r_stop, c_start, c_stop = bounds
    #
    x = g.x[r_start:r_stop, c_start:c_stop]
    y = g.y[r_start:r_stop, c_start:c_stop]
    #
    h_n = g.h_n[r_start:r_stop, c_start:c_stop].copy()
    B_n = g.B_n[r_start:r_stop, c_start:c_stop].copy()
    #
    te_n = g.t_erupted[r_start:r_stop, c_start:c_stop].copy()
    ti_n = g.t_inundation[r_start:r_stop, c_start:c_stop].copy()
    abs_U_s = g.abs_Usurf[r_start:r_stop, c_start:c_stop].copy()
    #
    return x, y, h_n, B_n, te_n, ti_n, abs_U_s

def local_mean(grid):
    out = grid.copy()
    C = out[1:-1,1:-1]
    W = out[1:-1,0:-2]
    E = out[1:-1,2:]
    N = out[0:-2,1:-1]
    S = out[2:,1:-1]
    NW = out[0:-2,0:-2]
    SE = out[2:,2:]
    NE = out[2:,0:-2]
    SW = out[0:-2,2:]
    #out[1:-1,1:-1] = (C + W + E + N + S + NW + SE + NE + SW) / 9.
    out[1:-1,1:-1] = ne.evaluate('C + W + E + N + S + NW + SE + NE + SW') / 9.
    return out

def local_max(grid):
    out = grid.copy()
    C = out[1:-1,1:-1]
    W = out[1:-1,0:-2]
    E = out[1:-1,2:]
    N = out[0:-2,1:-1]
    S = out[2:,1:-1]
    NW = out[0:-2,0:-2]
    SE = out[2:,2:]
    NE = out[2:,0:-2]
    SW = out[0:-2,2:]
    out[1:-1,1:-1] = np.maximum.reduce([C, W, E, N, S, NW, SE, NE, SW])
    return out

def local_and(grid):
    # grid is logical array
    # output: true if 9-cell neighborhood all true
    out = grid.copy()
    C = out[1:-1,1:-1]
    W = out[1:-1,0:-2]
    E = out[1:-1,2:]
    N = out[0:-2,1:-1]
    S = out[2:,1:-1]
    NW = out[0:-2,0:-2]
    SE = out[2:,2:]
    NE = out[2:,0:-2]
    SW = out[0:-2,2:]
    out[1:-1,1:-1] = np.logical_and.reduce([C, W, E, N, S, NW, SE, NE, SW])
    return out

def local_max_EWNS(grid):
    out = grid.copy()
    C = out[1:-1,1:-1]
    W = out[1:-1,0:-2]
    E = out[1:-1,2:]
    N = out[0:-2,1:-1]
    S = out[2:,1:-1]
    out[1:-1,1:-1] = np.maximum.reduce([C, W, E, N, S])
    return out

def local_and_EWNS(grid):
    # grid is logical array
    # output: true if neighborhood all true
    out = grid.copy()
    C = out[1:-1,1:-1]
    W = out[1:-1,0:-2]
    E = out[1:-1,2:]
    N = out[0:-2,1:-1]
    S = out[2:,1:-1]
    out[1:-1,1:-1] = np.logical_and.reduce([C, W, E, N, S])
    return out

def active_indicator(h, q, buffer):
    #I0 = np.logical_or(h > 0, q != 0)
    I0 = h + q != 0
    if I0.any():
        I_num = I0.astype(h.dtype)
        I_num = local_mean(I_num)
        return (I_num > 0.)
    else:
        I0[2:-2,2:-2] = True
        return I0

################################################################################
# Grid transfer
################################################################################
def interp(coarse):
    sh_coarse = np.array(coarse.shape)
    sh_fine = 2 * (sh_coarse-1) + 1
    fine = np.zeros(sh_fine, dtype = coarse.dtype)
    #
    fine[::2,::2] = coarse.copy() # injection
    fine[::2,1:-1:2] = 0.5 * (coarse[:,0:-1] + coarse[:,1:])  # sweep cols
    fine[1:-1:2,::2] = 0.5 * (coarse[0:-1,:] + coarse[1:,:])  # sweep rows
    fine[1:-1:2,1:-1:2] = 0.25 * (coarse[0:-1,0:-1] + coarse[1:,1:] + coarse[0:-1,1:] + coarse[1:,0:-1]) # fill holes
    #
    return fine

def restrict(fine):
    # Full Weighting "Molecule"
    coarse = fine[::2,::2].copy() # captures edges
    # Fill Interior
    coarse[1:-1,1:-1] = (
        1./4  *  fine[2:-2:2,2:-2:2] +
        1./8  * (fine[2:-2:2,1:-3:2] + fine[2:-2:2,3:-1:2] + fine[1:-3:2,2:-2:2] + fine[3:-1:2,2:-2:2]) +
        1./16 * (fine[1:-3:2,1:-3:2] + fine[3:-1:2,3:-1:2] + fine[3:-1:2,1:-3:2] + fine[1:-3:2,3:-1:2])
        )
    return coarse

################################################################################
# MUSCL reconstruction
################################################################################

def flux_lim_van_Leer(d1,d2):
    # d1 = u_{i} - u_{i-1}
    # d2 = u_{i+1} - u_{i}
    # r = d1 / d2
    #
    f = 2 + np.zeros(d1.shape)
    valid = abs(d2) > p.pos_eps
    r = d1[valid] / d2[valid]
    abs_r = abs(r)
    f[valid] = (r + abs_r) / (1 + abs_r)
    #
    return f


def flux_lim_superbee(d1,d2):
    # d1 = u_{i} - u_{i-1}
    # d2 = u_{i+1} - u_{i}
    # r = d1 / d2
    #
    f = 2 + np.zeros(d1.shape)
    valid = abs(d2) > p.pos_eps
    r = d1[valid] / d2[valid]
    f[valid] = np.maximum.reduce([np.zeros(r.shape), np.minimum(2*r, 1), np.minimum(r, 2)])
    #
    return f


def muscl_linear(u, I, flux_lim):
    # --------------------------------------------------------------------------
    # MUSCL Scheme: Piecewise-Linear Reconstruction
    #
    I_ij = I[2:-2,2:-2]
    u_NN = u[0:-4,2:-2][I_ij]
    u_N  = u[1:-3,2:-2][I_ij]
    u_W  = u[2:-2,1:-3][I_ij]
    u_WW = u[2:-2,0:-4][I_ij]
    u_ij = u[I] # Note: u[I] == u[2:-2,2:-2][I_ij]
    u_EE = u[2:-2,4:][I_ij]
    u_E  = u[2:-2,3:-1][I_ij]
    u_S  = u[3:-1,2:-2][I_ij]
    u_SS = u[4:,2:-2][I_ij]
    #
    du_WW = u_W - u_WW # (rows-4, cols-4)
    du_W  = u_ij - u_W  # (rows-4, cols-4)
    du_E  = u_E - u_ij  # (rows-4, cols-4)
    du_EE = u_EE - u_E  # (rows-4, cols-4)
    #
    du_NN = u_NN - u_N  # (rows-4, cols-4)
    du_N  = u_N - u_ij  # (rows-4, cols-4)
    du_S  = u_ij - u_S  # (rows-4, cols-4)
    du_SS = u_S - u_SS  # (rows-4, cols-4)
    #
    half_lim_WE_du_E = 0.5 * flux_lim(du_W, du_E) * du_E
    half_lim_SN_du_N = 0.5 * flux_lim(du_S, du_N) * du_N
    #
    u_E_R = u_E  - 0.5 * flux_lim(du_E,du_EE) * du_EE # East face, Right state
    u_E_L = u_ij + half_lim_WE_du_E                            # East face, Left state
    u_W_R = u_ij - half_lim_WE_du_E                            # West face, Right state
    u_W_L = u_W  + 0.5 * flux_lim(du_WW,du_W) * du_W  # West face, Left state
    #
    u_N_U = u_N  - 0.5 * flux_lim(du_N,du_NN) * du_NN # North face, Upper state
    u_N_D = u_ij + half_lim_SN_du_N                            # North face, Lower state
    u_S_U = u_ij - half_lim_SN_du_N                            # South face, Upper state
    u_S_D = u_S  + 0.5 * flux_lim(du_SS,du_S) * du_S  # South face, Lower state
    #
    # --------------------------------------------------------------------------
    #
    return u_E_R, u_E_L, u_W_R, u_W_L, u_N_U, u_N_D, u_S_U, u_S_D

################################################################################
# timesteps
################################################################################

def dt_muscl(x_speed_max, y_speed_max, K_max, h_max, q_max):
    if (h_max <= p.tiny_flow):
        dt = p.tiny_flow / (q_max + p.pos_eps) + 1e-3
        dt_type = 'init'
    else:
        rate_hyperbolic = max(x_speed_max/p.dx, y_speed_max/p.dy)
        rate_parabolic = (1./p.dx**2 + 1./p.dy**2) * K_max
        dt_hyperbolic  = p.cfl_max / (rate_hyperbolic + p.pos_eps)
        dt_parabolic  = p.lambda_max / (rate_parabolic + p.pos_eps)
        #
        if dt_hyperbolic < dt_parabolic:
            dt = p.safety_factor * dt_hyperbolic
            dt_type = 'hyp'
        else:
            dt = p.safety_factor * dt_parabolic
            dt_type = 'par'
            #
        #
    #
    if dt > p.dt_max:
        dt = p.dt_max
        dt_type = 'max'
    #
    return dt, dt_type

def dt_hyperbolic(x_speed_max, y_speed_max, h_max, q_max):
    if (h_max <= p.tiny_flow):
        dt = p.tiny_flow / (q_max + p.pos_eps) + 1e-3
        dt_type = 'init'
    else:
        rate_hyperbolic = max(x_speed_max/p.dx, y_speed_max/p.dy)
        dt = p.safety_factor * p.cfl_max / (rate_hyperbolic + p.pos_eps)
        dt_type = 'hyp'
        #
    #
    if dt > p.dt_max:
        dt = p.dt_max
        dt_type = 'max'
    #
    return dt, dt_type


################################################################################
# KT 2000 MUSCL Step
################################################################################
def step_KT2000_FD2_RFA(h, B_n, te, ti, q_n, abs_U_s, t_n, I):
    # --------------------------------------------------------------------------
    # methods:
    # physics_quality   = 'research' (R)
    # h_update          = 'full' (F)
    # te_update         = 'advection' (A)
    # --------------------------------------------------------------------------
    # Compute all fields here:
    phi_S = therm.surface_BL(t_n, te, h, p.core_temperature)
    cryst_core = rheo.cryst_avrami(t_n, te, I)
    R_n, Rsurf_n, fluidity_n, abs_grad_S_n, jeffreys_efficiency = rheo.rheo_factor_bl_S_all_outputs(h, B_n, phi_S, cryst_core)
    dBdt = rheo.dBdt(t_n, ti, h, p.core_temperature, abs_U_s, fluidity_n, q_n)
    #
    # --------------------------------------------------------------------------
    # Freezout:
    topo.freeze(h, B_n, dBdt, ti, phi_S, cryst_core, jeffreys_efficiency, t_n)
    #
    # --------------------------------------------------------------------------
    # Active subsetting # Note: Q[I] == Q[2:-2,2:-2][I_ij]
    K_n = ne.evaluate('R_n * h**3') # much faster cubing
    I_ij = I[2:-2,2:-2] # (rows-4, cols-4)
    q_ij = q_n[I] # (n_active,)
    h_ij = h[I] # (n_active,)
    B_ij = B_n[I] # (n_active,)
    K_ij = K_n[I] # (n_active,)
    R_ij = R_n[I] # (n_active,)
    h_RLUD = np.array([h[2:-2,3:-1][I_ij], h[2:-2,1:-3][I_ij], h[1:-3,2:-2][I_ij], h[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    B_RLUD = np.array([B_n[2:-2,3:-1][I_ij], B_n[2:-2,1:-3][I_ij], B_n[1:-3,2:-2][I_ij], B_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    K_RLUD = np.array([K_n[2:-2,3:-1][I_ij], K_n[2:-2,1:-3][I_ij], K_n[1:-3,2:-2][I_ij], K_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    R_RLUD = np.array([R_n[2:-2,3:-1][I_ij], R_n[2:-2,1:-3][I_ij], R_n[1:-3,2:-2][I_ij], R_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    #
    # --------------------------------------------------------------------------
    # Edge Variables
    K = 0.5 * (K_ij[None,:] + K_RLUD) # (4, n_active) # order: E, W, N, S
    R = 0.5 * (R_ij[None,:] + R_RLUD) # (4, n_active) # order: E, W, N, S
    dxy = np.array([p.dx, -p.dx, p.dy, -p.dy])
    grad_h = (h_RLUD - h_ij[None,:]) / dxy[:,None] # (4, n_active) # order: E, W, N, S
    neg_grad_B = (B_ij[None,:] - B_RLUD) / dxy[:,None] # (4, n_active) # order: E, W, N, S
    neg_R_grad_B = R * neg_grad_B
    #
    # --------------------------------------------------------------------------
    # MUSCL reconstruction of edge states (piecewise linear)
    all_h_states = np.array(muscl_linear(h, I,flux_lim_van_Leer)) # (8, n_active)
    h_RRUU = all_h_states[0::2] # (4, n_active) RU: (EWNS, n_active)
    h_LLDD = all_h_states[1::2] # (4, n_active) LD: (EWNS, n_active)
    h2_RRUU = h_RRUU**2 # (4, n_active) RU: (EWNS, n_active)
    h2_LLDD = h_LLDD**2 # (4, n_active) LD: (EWNS, n_active)
    #
    # --------------------------------------------------------------------------
    # CFL related speeds
    # a: dF/dh ~ 3 * R * -grad(B) * h**2 == 3 * Actual Advection Speed
    max_h2 = np.maximum(h2_RRUU, h2_LLDD) # (4, n_active) # order: E, W, N, S
    a = abs(neg_R_grad_B) * 3*max_h2 # (4, n_active) # order: E, W, N, S
    amax = np.max(a, axis = -1) # (4,) # order: E, W, N, S
    x_speed_max = amax[0:2].max()
    y_speed_max = amax[2:4].max()
    #
    # --------------------------------------------------------------------------
    # Time step Calculation:
    dt, dt_type = dt_muscl(x_speed_max, y_speed_max, K_ij.max(), h_ij.max(), q_ij.max())
    #
    # --------------------------------------------------------------------------
    # Fluxes from Kurganov and Tadmor (2000)
    all_fluxes_h = 0.5 * (
        (neg_R_grad_B*h2_RRUU - a)*h_RRUU +
        (neg_R_grad_B*h2_LLDD + a)*h_LLDD
        ) - K * grad_h # (4, n_active)
    #
    # --------------------------------------------------------------------------
    # Logic to stop uphill fluxes
    neg_grad_S = neg_grad_B - grad_h # order: E, W, N, S
    is_uphill_flow = np.sign(neg_grad_S * all_fluxes_h) == -1 # (4, n_active) # order: E, W, N, S
    all_fluxes_h[is_uphill_flow] = 0 # (4, n_active) # order: E, W, N, S
    #
    # --------------------------------------------------------------------------
    # Forward Euler Update
    h[I] += dt * (q_ij - np.sum(all_fluxes_h / dxy[:,None], axis = 0))
    substeps = 1
    #
    # --------------------------------------------------------------------------
    # update t_erupted, get surface speed
    abs_Usurf = update_te_from_advection(te, Rsurf_n, neg_grad_S, h2_RRUU, h2_LLDD, max_h2, q_n, dxy, t_n, dt, is_uphill_flow, I, I_ij) # Implicit globals: Rsurf_n, neg_grad_S, h2_RRUU, h2_LLDD, max_h2, q_n, dxy, t_n, dt, is_uphill_flow, I, I_ij
    # --------------------------------------------------------------------------
    # Update Base
    dB = np.minimum(dt * dBdt, h) # at most can freeze h
    h -= dB
    B_n += dB
    #
    # --------------------------------------------------------------------------
    return h, te, abs_Usurf, B_n, dt, dt_type, substeps



def step_KT2000_FD2_OHA(h, B_n, te, ti, q_n, abs_U_s, t_n, I):
    # --------------------------------------------------------------------------
    # methods:
    # physics_quality   = 'operations' (O)
    # h_update          = 'hyperbolic_approx' (H)
    # te_update         = 'advection' (A)
    # --------------------------------------------------------------------------
    # Compute all fields here:
    phi_S = therm.surface_BL(t_n, te, h, p.core_temperature)
    #phi_S = therm.surface_BL_simple(t_n, te, h, p.core_temperature)
    cryst_core = rheo.cryst_avrami(t_n, te, I)
    R_n, Rsurf_n, fluidity_n, abs_grad_S_n, jeffreys_efficiency = rheo.rheo_factor_bl_S_all_outputs(h, B_n, phi_S, cryst_core)
    dBdt = rheo.dBdt(t_n, ti, h, p.core_temperature, abs_U_s, fluidity_n, q_n)
    #
    # --------------------------------------------------------------------------
    # Freezout:
    topo.freeze(h, B_n, dBdt, ti, phi_S, cryst_core, jeffreys_efficiency, t_n)
    #
    # --------------------------------------------------------------------------
    # Active subsetting # Note: Q[I] == Q[2:-2,2:-2][I_ij]
    S_n = B_n + h
    K_n = ne.evaluate('R_n * h**3') # much faster cubing
    I_ij = I[2:-2,2:-2] # (rows-4, cols-4)
    q_ij = q_n[I] # (n_active,)
    #
    R_ij = R_n[I] # (n_active,)
    S_ij = S_n[I] # (n_active,)
    R_RLUD = np.array([R_n[2:-2,3:-1][I_ij], R_n[2:-2,1:-3][I_ij], R_n[1:-3,2:-2][I_ij], R_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    S_RLUD = np.array([S_n[2:-2,3:-1][I_ij], S_n[2:-2,1:-3][I_ij], S_n[1:-3,2:-2][I_ij], S_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    # --------------------------------------------------------------------------
    # Edge Variables
    R = 0.5 * (R_ij[None,:] + R_RLUD) # (4, n_active) # order: E, W, N, S
    dxy = np.array([p.dx, -p.dx, p.dy, -p.dy]) # (4,) # order: E, W, N, S
    neg_grad_S = (S_ij[None,:] - S_RLUD) / dxy[:,None] # (4, n_active) # order: E, W, N, S
    neg_R_grad_S = R * neg_grad_S
    #
    # --------------------------------------------------------------------------
    # MUSCL reconstruction of edge states (piecewise linear)
    all_h_states = np.array(muscl_linear(h, I, flux_lim_van_Leer)) # (8, n_active)
    h_RRUU = all_h_states[0::2] # (4, n_active) RU: (EWNS, n_active)
    h_LLDD = all_h_states[1::2] # (4, n_active) LD: (EWNS, n_active)
    h2_RRUU = h_RRUU**2 # (4, n_active) RU: (EWNS, n_active)
    h2_LLDD = h_LLDD**2 # (4, n_active) LD: (EWNS, n_active)
    #
    # --------------------------------------------------------------------------
    # CFL related speeds
    # dF/du = 3 * R * -grad(S) * h**2 == 3 * Actual Advection Speed
    max_h2 = np.maximum(h2_RRUU, h2_LLDD) # (4, n_active) # order: E, W, N, S
    a = abs(neg_R_grad_S) * 3*max_h2 # (4, n_active) # order: E, W, N, S
    amax = np.max(a, axis = -1) # (4,) # order: E, W, N, S
    x_speed_max = amax[0:2].max()
    y_speed_max = amax[2:4].max()
    #
    Kmax = np.max(K_n[I])
    #
    # --------------------------------------------------------------------------
    # Timestep
    # set dt: either hyperbolic, initial, or max_val
    dt, dt_type = dt_muscl(x_speed_max, y_speed_max, Kmax, h.max(), q_ij.max())
    #
    # --------------------------------------------------------------------------
    # Fluxes from Kurganov and Tadmor (2000)
    all_convective_fluxes = 0.5 * (
        (neg_R_grad_S*h2_RRUU - a)*h_RRUU +
        (neg_R_grad_S*h2_LLDD + a)*h_LLDD
        ) # (4, n_active) # convective fluxes
    #
    # ----------------------------------------------------------------------
    # Logic to stop uphill fluxes
    is_uphill_flow = np.sign(neg_grad_S * all_convective_fluxes) == -1 # (4, n_active) # order: E, W, N, S
    all_convective_fluxes[is_uphill_flow] = 0 # (4, n_active) # order: E, W, N, S
    #
    # ----------------------------------------------------------------------
    # Forward Euler Update
    h[I] += dt * (q_ij - np.sum(all_convective_fluxes/dxy[:,None], axis = 0))
    substeps = 1
    # --------------------------------------------------------------------------
    # update t_erupted, get surface speed
    abs_Usurf = update_te_from_advection(te, Rsurf_n, neg_grad_S, h2_RRUU, h2_LLDD, max_h2, q_n, dxy, t_n, dt, is_uphill_flow, I, I_ij) # Implicit globals: Rsurf_n, neg_grad_S, h2_RRUU, h2_LLDD, max_h2, q_n, dxy, t_n, dt, is_uphill_flow, I, I_ij
    #
    # --------------------------------------------------------------------------
    # Update Base
    dB = np.minimum(dt * dBdt, h) # at most can freeze h
    h -= dB
    B_n += dB
    #
    return h, te, abs_Usurf, B_n, dt, dt_type, substeps



def step_KT2000_FD2_OSA(h, B_n, te, ti, q_n, abs_U_s, t_n, I):
    # --------------------------------------------------------------------------
    # methods:
    # physics_quality   = 'operations' (O)
    # h_update          = 'parabolic_substeps' (S)
    # te_update         = 'advection' (A)
    # --------------------------------------------------------------------------
    # Compute all fields here:
    phi_S = therm.surface_BL_simple(t_n, te, h, p.core_temperature)
    cryst_core = rheo.cryst_avrami(t_n, te, I)
    R_n, Rsurf_n, fluidity_n, abs_grad_S_n, jeffreys_efficiency = rheo.rheo_factor_bl_S_all_outputs(h, B_n, phi_S, cryst_core)
    dBdt = rheo.dBdt(t_n, ti, h, p.core_temperature, abs_U_s, fluidity_n, q_n)
    #
    # --------------------------------------------------------------------------
    # Freezout:
    topo.freeze(h, B_n, dBdt, ti, phi_S, cryst_core, jeffreys_efficiency, t_n)
    #
    # --------------------------------------------------------------------------
    # Active subsetting # Note: Q[I] == Q[2:-2,2:-2][I_ij]
    K_n = ne.evaluate('R_n * h**3') # much faster cubing
    I_ij = I[2:-2,2:-2] # (rows-4, cols-4)
    q_ij = q_n[I] # (n_active,)
    #
    R_ij = R_n[I] # (n_active,)
    K_ij = K_n[I] # (n_active,)
    B_ij = B_n[I] # (n_active,)
    R_RLUD = np.array([R_n[2:-2,3:-1][I_ij], R_n[2:-2,1:-3][I_ij], R_n[1:-3,2:-2][I_ij], R_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    K_RLUD = np.array([K_n[2:-2,3:-1][I_ij], K_n[2:-2,1:-3][I_ij], K_n[1:-3,2:-2][I_ij], K_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    B_RLUD = np.array([B_n[2:-2,3:-1][I_ij], B_n[2:-2,1:-3][I_ij], B_n[1:-3,2:-2][I_ij], B_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    #
    # --------------------------------------------------------------------------
    # Edge Variables
    R = 0.5 * (R_ij[None,:] + R_RLUD) # (4, n_active) # order: E, W, N, S
    K = 0.5 * (K_ij[None,:] + K_RLUD) # (4, n_active) # order: E, W, N, S
    dxy = np.array([p.dx, -p.dx, p.dy, -p.dy]) # (4,) # order: E, W, N, S
    neg_grad_B = (B_ij[None,:] - B_RLUD) / dxy[:,None] # (4, n_active) # order: E, W, N, S
    neg_R_grad_B = R * neg_grad_B
    #
    # --------------------------------------------------------------------------
    # MUSCL reconstruction of edge states (piecewise linear)
    all_h_states = np.array(muscl_linear(h, I, flux_lim_van_Leer)) # (8, n_active)
    h_RRUU = all_h_states[0::2] # (4, n_active) RU: (EWNS, n_active)
    h_LLDD = all_h_states[1::2] # (4, n_active) LD: (EWNS, n_active)
    h2_RRUU = h_RRUU**2 # (4, n_active) RU: (EWNS, n_active)
    h2_LLDD = h_LLDD**2 # (4, n_active) LD: (EWNS, n_active)
    #
    # --------------------------------------------------------------------------
    # CFL related speeds
    # dF/du = 3 * R * -grad(B) * h**2 == 3 * Actual Advection Speed
    max_h2 = np.maximum(h2_RRUU, h2_LLDD) # (4, n_active) # order: E, W, N, S
    a = abs(neg_R_grad_B) * 3*max_h2 # (4, n_active) # order: E, W, N, S
    amax = np.max(a, axis = -1) # (4,) # order: E, W, N, S
    x_speed_max = amax[0:2].max()
    y_speed_max = amax[2:4].max()
    #
    # --------------------------------------------------------------------------
    # Fluxes from Kurganov and Tadmor (2000)
    all_convective_fluxes_h = 0.5 * (
        (neg_R_grad_B*h2_RRUU - a)*h_RRUU +
        (neg_R_grad_B*h2_LLDD + a)*h_LLDD
        ) # (4, n_active)
    #
    # --------------------------------------------------------------------------
    # Timestep
    # set dt: either hyperbolic, initial, or max_val
    dt, dt_type = dt_hyperbolic(x_speed_max, y_speed_max, h.max(), q_ij.max())
    dt_ftcs = p.safety_factor * p.lambda_max / ((K.max() + p.pos_eps) * (1/p.dx**2 + 1/p.dy**2))
    #
    # --------------------------------------------------------------------------
    # Iterate on substeps for diffusion here (using FTCS):
    # all_convective_fluxes, K, q stay constant while Diffusion is carried forward
    h_sub = h.copy()
    substeps = 0
    tau_sub = 0 # measures progress through dt interval: t_n <= t_n + tau_sub <= t_n + dt
    while tau_sub < dt:
        # ----------------------------------------------------------------------
        # diffusion-limited, convection-limited, or time-remaining in dt
        dt_sub = min(dt_ftcs, dt, dt-tau_sub)
        #
        # ----------------------------------------------------------------------
        # Compute grad(h), total fluxes
        h_ij = h_sub[I] # (n_active,)
        h_RLUD = np.array([h_sub[2:-2,3:-1][I_ij], h_sub[2:-2,1:-3][I_ij], h_sub[1:-3,2:-2][I_ij], h_sub[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
        grad_h = (h_RLUD - h_ij[None,:]) / dxy[:,None] # (4, n_active) # order: E, W, N, S
        all_fluxes_h = all_convective_fluxes_h - K*grad_h # (4, n_active) # order: E, W, N, S
        #
        # ----------------------------------------------------------------------
        # Logic to stop uphill fluxes
        neg_grad_S = neg_grad_B - grad_h # order: E, W, N, S
        is_uphill_flow = np.sign(neg_grad_S * all_fluxes_h) == -1 # (4, n_active) # order: E, W, N, S
        all_fluxes_h[is_uphill_flow] = 0 # (4, n_active) # order: E, W, N, S
        #
        # ----------------------------------------------------------------------
        # Forward Euler Update
        h_sub[I] += dt_sub * (q_ij - np.sum(all_fluxes_h/dxy[:,None], axis = 0))
        tau_sub += dt_sub
        substeps += 1
        #
        # ----------------------------------------------------------------------
    # h_sub is now approximately h_{n+1}
    h = h_sub.copy()
    #
    # --------------------------------------------------------------------------
    # update t_erupted, get surface speed
    abs_Usurf = update_te_from_advection(te, Rsurf_n, neg_grad_S, h2_RRUU, h2_LLDD, max_h2, q_n, dxy, t_n, dt, is_uphill_flow, I, I_ij)
    #
    # --------------------------------------------------------------------------
    # Update Base
    dB = np.minimum(dt * dBdt, h) # at most can freeze h
    h -= dB
    B_n += dB
    #
    # --------------------------------------------------------------------------
    return h, te, abs_Usurf, B_n, dt, dt_type, substeps



def step_KT2000_FD2_OHI(h, B_n, te, ti, q_n, abs_U_s, t_n, I):
    # --------------------------------------------------------------------------
    # methods:
    # physics_quality   = 'operations' (O)
    # h_update          = 'hyperbolic_approx' (H)
    # te_update         = 'from_inundation' (I)
    # --------------------------------------------------------------------------
    # Compute all fields here:
    phi_S = therm.surface_BL_simple(t_n, te, h, p.core_temperature)
    cryst_core = rheo.cryst_avrami(t_n, te, I)
    R_n, Rsurf_n, fluidity_n, abs_grad_S_n, jeffreys_efficiency = rheo.rheo_factor_bl_S_all_outputs(h, B_n, phi_S, cryst_core)
    dBdt = rheo.dBdt(t_n, ti, h, p.core_temperature, abs_U_s, fluidity_n, q_n)
    #
    # --------------------------------------------------------------------------
    # update t_erupted, get surface speed
    abs_Usurf = update_te_from_inundation(te, Rsurf_n, abs_grad_S_n, h, R_n, ti, t_n, dt, I) # Implicit Globals: Rsurf_n, abs_grad_S_n, h, R_n, ti, t_n, dt, I
    #
    # --------------------------------------------------------------------------
    # Freezout:
    topo.freeze(h, B_n, dBdt, ti, phi_S, cryst_core, jeffreys_efficiency, t_n)
    #
    # --------------------------------------------------------------------------
    # Active subsetting # Note: Q[I] == Q[2:-2,2:-2][I_ij]
    S_n = B_n + h
    K_n = ne.evaluate('R_n * h**3') # much faster cubing
    I_ij = I[2:-2,2:-2] # (rows-4, cols-4)
    q_ij = q_n[I] # (n_active,)
    #
    R_ij = R_n[I] # (n_active,)
    S_ij = S_n[I] # (n_active,)
    R_RLUD = np.array([R_n[2:-2,3:-1][I_ij], R_n[2:-2,1:-3][I_ij], R_n[1:-3,2:-2][I_ij], R_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    S_RLUD = np.array([S_n[2:-2,3:-1][I_ij], S_n[2:-2,1:-3][I_ij], S_n[1:-3,2:-2][I_ij], S_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    # --------------------------------------------------------------------------
    # Edge Variables
    R = 0.5 * (R_ij[None,:] + R_RLUD) # (4, n_active) # order: E, W, N, S
    dxy = np.array([p.dx, -p.dx, p.dy, -p.dy]) # (4,) # order: E, W, N, S
    neg_grad_S = (S_ij[None,:] - S_RLUD) / dxy[:,None] # (4, n_active) # order: E, W, N, S
    neg_R_grad_S = R * neg_grad_S
    #
    # --------------------------------------------------------------------------
    # MUSCL reconstruction of edge states (piecewise linear)
    all_h_states = np.array(muscl_linear(h, I, flux_lim_van_Leer)) # (8, n_active)
    h_RRUU = all_h_states[0::2] # (4, n_active) RU: (EWNS, n_active)
    h_LLDD = all_h_states[1::2] # (4, n_active) LD: (EWNS, n_active)
    h2_RRUU = h_RRUU**2 # (4, n_active) RU: (EWNS, n_active)
    h2_LLDD = h_LLDD**2 # (4, n_active) LD: (EWNS, n_active)
    #
    # --------------------------------------------------------------------------
    # CFL related speeds
    # dF/du = 3 * R * -grad(S) * h**2 == 3 * Actual Advection Speed
    a = abs(neg_R_grad_S) * 3*np.maximum(h2_RRUU, h2_LLDD) # (4, n_active) # order: E, W, N, S
    amax = np.max(a, axis = -1) # (4,) # order: E, W, N, S
    x_speed_max = amax[0:2].max()
    y_speed_max = amax[2:4].max()
    #
    # --------------------------------------------------------------------------
    # Timestep
    # set dt: either hyperbolic, initial, or max_val
    dt, dt_type = dt_muscl(x_speed_max, y_speed_max, K_n.max(), h.max(), q_ij.max())
    #
    # --------------------------------------------------------------------------
    # Fluxes from Kurganov and Tadmor (2000)
    all_convective_fluxes = 0.5 * (
        (neg_R_grad_S*h2_RRUU - a)*h_RRUU +
        (neg_R_grad_S*h2_LLDD + a)*h_LLDD
        ) # (4, n_active) # convective fluxes
    #
    # ----------------------------------------------------------------------
    # Logic to stop uphill fluxes
    is_uphill_flow = np.sign(neg_grad_S * all_convective_fluxes) == -1 # (4, n_active) # order: E, W, N, S
    all_convective_fluxes[is_uphill_flow] = 0 # (4, n_active) # order: E, W, N, S
    #
    # ----------------------------------------------------------------------
    # Forward Euler Update
    h[I] += dt * (q_ij - np.sum(all_convective_fluxes/dxy[:,None], axis = 0))
    substeps = 1
    # --------------------------------------------------------------------------
    # Update Base
    dB = np.minimum(dt * dBdt, h) # at most can freeze h
    h -= dB
    B_n += dB
    #
    return h, te, abs_Usurf, B_n, dt, dt_type, substeps



def step_KT2000_FD2_OSI(h, B_n, te, ti, q_n, abs_U_s, t_n, I):
    # --------------------------------------------------------------------------
    # methods:
    # physics_quality   = 'operations' (O)
    # h_update          = 'parabolic_substeps' (S)
    # te_update         = 'from_inundation' (I)
    # --------------------------------------------------------------------------
    # Compute all fields here:
    phi_S = therm.surface_BL_simple(t_n, te, h, p.core_temperature)
    cryst_core = rheo.cryst_avrami(t_n, te, I)
    R_n, Rsurf_n, fluidity_n, abs_grad_S_n, jeffreys_efficiency = rheo.rheo_factor_bl_S_all_outputs(h, B_n, phi_S, cryst_core)
    dBdt = rheo.dBdt(t_n, ti, h, p.core_temperature, abs_U_s, fluidity_n, q_n)
    #
    # --------------------------------------------------------------------------
    # update t_erupted, get surface speed
    abs_Usurf = update_te_from_inundation(te, Rsurf_n, abs_grad_S_n, h, R_n, ti, t_n, dt, I) # Implicit Globals: Rsurf_n, abs_grad_S_n, h, R_n, ti, t_n, dt, I
    #
    # --------------------------------------------------------------------------
    # Freezout:
    topo.freeze(h, B_n, dBdt, ti, phi_S, cryst_core, jeffreys_efficiency, t_n)
    #
    # --------------------------------------------------------------------------
    # Active subsetting # Note: Q[I] == Q[2:-2,2:-2][I_ij]
    K_n = ne.evaluate('R_n * h**3') # much faster cubing
    I_ij = I[2:-2,2:-2] # (rows-4, cols-4)
    q_ij = q_n[I] # (n_active,)
    #
    R_ij = R_n[I] # (n_active,)
    K_ij = K_n[I] # (n_active,)
    B_ij = B_n[I] # (n_active,)
    R_RLUD = np.array([R_n[2:-2,3:-1][I_ij], R_n[2:-2,1:-3][I_ij], R_n[1:-3,2:-2][I_ij], R_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    K_RLUD = np.array([K_n[2:-2,3:-1][I_ij], K_n[2:-2,1:-3][I_ij], K_n[1:-3,2:-2][I_ij], K_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    B_RLUD = np.array([B_n[2:-2,3:-1][I_ij], B_n[2:-2,1:-3][I_ij], B_n[1:-3,2:-2][I_ij], B_n[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
    # --------------------------------------------------------------------------
    # Edge Variables
    R = 0.5 * (R_ij[None,:] + R_RLUD) # (4, n_active) # order: E, W, N, S
    K = 0.5 * (K_ij[None,:] + K_RLUD) # (4, n_active) # order: E, W, N, S
    dxy = np.array([p.dx, -p.dx, p.dy, -p.dy]) # (4,) # order: E, W, N, S
    neg_grad_B = (B_ij[None,:] - B_RLUD) / dxy[:,None] # (4, n_active) # order: E, W, N, S
    neg_R_grad_B = R * neg_grad_B
    #
    # --------------------------------------------------------------------------
    # MUSCL reconstruction of edge states (piecewise linear)
    all_h_states = np.array(muscl_linear(h, I, flux_lim_van_Leer)) # (8, n_active)
    h_RRUU = all_h_states[0::2] # (4, n_active) RU: (EWNS, n_active)
    h_LLDD = all_h_states[1::2] # (4, n_active) LD: (EWNS, n_active)
    h2_RRUU = h_RRUU**2 # (4, n_active) RU: (EWNS, n_active)
    h2_LLDD = h_LLDD**2 # (4, n_active) LD: (EWNS, n_active)
    #
    # --------------------------------------------------------------------------
    # CFL related speeds
    # dF/du = 3 * R * -grad(B) * h**2 == 3 * Actual Advection Speed
    a = abs(neg_R_grad_B) * 3*np.maximum(h2_RRUU, h2_LLDD) # (4, n_active) # order: E, W, N, S
    amax = np.max(a, axis = -1) # (4,) # order: E, W, N, S
    x_speed_max = amax[0:2].max()
    y_speed_max = amax[2:4].max()
    #
    # --------------------------------------------------------------------------
    # Fluxes from Kurganov and Tadmor (2000)
    all_convective_fluxes = 0.5 * (
        (neg_R_grad_B*h2_RRUU - a)*h_RRUU +
        (neg_R_grad_B*h2_LLDD + a)*h_LLDD
        ) # (4, n_active) # convective fluxes
    #
    # --------------------------------------------------------------------------
    # Timestep
    # set dt: either hyperbolic, initial, or max_val
    dt, dt_type = dt_hyperbolic(x_speed_max, y_speed_max, h.max(), q_ij.max())
    dt_ftcs = p.safety_factor * p.lambda_max / ((K.max() + p.pos_eps) * (1/p.dx**2 + 1/p.dy**2))
    #
    # --------------------------------------------------------------------------
    # Iterate on substeps for diffusion here (using FTCS):
    # all_convective_fluxes, K, q stay constant while Diffusion is carried forward
    h_sub = h.copy()
    substeps = 0
    tau_sub = 0 # measures progress through dt interval: t_n <= t_n + tau_sub <= t_n + dt
    while tau_sub < dt:
        # ----------------------------------------------------------------------
        # diffusion-limited, convection-limited, or time-remaining in dt
        dt_sub = min(dt_ftcs, dt, dt-tau_sub)
        #
        # ----------------------------------------------------------------------
        # Compute grad(h), total fluxes
        h_ij = h_sub[I] # (n_active,)
        h_RLUD = np.array([h_sub[2:-2,3:-1][I_ij], h_sub[2:-2,1:-3][I_ij], h_sub[1:-3,2:-2][I_ij], h_sub[3:-1,2:-2][I_ij]]) # (4, n_active) # order: R, W, N, S
        grad_h = (h_RLUD - h_ij[None,:]) / dxy[:,None] # (4, n_active) # order: E, W, N, S
        all_fluxes = all_convective_fluxes - K*grad_h # (4, n_active) # order: E, W, N, S
        #
        # ----------------------------------------------------------------------
        # Logic to stop uphill fluxes
        neg_grad_S = neg_grad_B - grad_h # order: E, W, N, S
        is_uphill_flow = np.sign(neg_grad_S * all_fluxes) == -1 # (4, n_active) # order: E, W, N, S
        all_fluxes[is_uphill_flow] = 0 # (4, n_active) # order: E, W, N, S
        #
        # ----------------------------------------------------------------------
        # Forward Euler Update
        h_sub[I] += dt_sub * (q_ij - np.sum(all_fluxes/dxy[:,None], axis = 0))
        tau_sub += dt_sub
        substeps += 1
        #
        # ----------------------------------------------------------------------
    # h_sub is now approximately h_{n+1}
    h = h_sub.copy()
    # alternative: instead of FTCS - use an implicit method and Jacobi's methods
    # difficulty: need to preserve logic to stop uphill fluxes
    #
    # --------------------------------------------------------------------------
    # Update t_erupted
    speed_ratio_ij = local_mean((R_n + p.pos_eps) / (Rsurf_n + p.pos_eps))[I]
    te_ij_update = t_n -  speed_ratio_ij * ti[I] # assume surface travel time is ~ 2/3 of bulk current
    te_ij_update[q_ij > 0] = t_n + dt # ensure source conditon  # (rows, cols)
    te[I] = te_ij_update
    #
    # --------------------------------------------------------------------------
    # Update Base
    dB = np.minimum(dt * dBdt, h) # at most can freeze h
    h -= dB
    B_n += dB
    #
    return h, te, abs_Usurf, B_n, dt, dt_type, substeps





################################################################################
#
################################################################################

import numpy as np
import numexpr as ne
import thermal as therm
from globals import params as p
################################################################################
# Rheological Model
################################################################################

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Specify Melt Viscosity Model
# must specify: viscosity(T), d_log_visc_dT(T), T_stiffened(T_core)

def melt_viscosity(T):
    # FLOWGO; matching VFT models <--> Dragoni (1989)
    #a = 0.04 # K-1
    #return  p.viscosity_vent * np.exp(a*(p.T_vent - T))
    #gv_contrast = np.log(1e12 / p.viscosity_vent)
    #dT = gv_contrast * (p.T_vent - p.T_glass) / (gv_contrast - a * (p.T_vent - p.T_glass))
    #return p.viscosity_vent * np.exp(-a*dT*(1 - dT/np.maximum(T - p.T_vent + dT, 0.1*dT)))
    # Lev et al 2012: visc(T) = 10**(-5.94 + 5500. / (T - 610.))
    a = d_log_visc_dT(p.core_temperature)
    return p.viscosity_vent * np.exp(a*(p.T_vent - T))
    #T0 = p.T_vent
    #mu0 = p.viscosity_vent
    #return ne.evaluate('mu0 * exp(a*(T0 - T))')

def d_log_visc_dT(T):
    # a =  -d(log(melt_visc))/dT
    # Lev et al 2012: visc(T) = 10**(-5.94 + 5500. / (T - 610.))
    # a_Lev = np.log(10) * 5500. * (T - 610.)**-2
    # a_lev_1050C = 0.0249
    # a_Dragoni = .04
    return np.log(10) * 5500. * (T - 610.)**-2

def T_stiffened(T_core):
    # Temp at which visc = 2*visc_core
    # assume all stiffening happens due to cooling of melt phase (no excess xtals)
    # ie: Temp at which:  visc_melt_stiff = 2*visc_melt_core
    # T_stiffened_Lev = 610. + 1 / (np.log10(2) / 5500. + 1 / (T_core - 610.))
    return T_core - np.log(2) / d_log_visc_dT(T_core)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def viscosity(T, phi):
    return np.maximum(1 - phi/p.phi_max, 1e-3)**-2.5 * melt_viscosity(T)

def fluidity(T, phi):
    return np.maximum(1 - phi/p.phi_max, 0)**2.5 / melt_viscosity(T)

def cryst_avrami(t, te, I):
    n = 4
    t_active = 1.0 * (t - te[I]) # (n_active,) # core arrival in terms of surface transit time
    C_vent = np.log(1. / (1. - p.cryst_vent))**(1./n)
    lam = (np.exp(n-1)/(n*(n-1)**(n-1)))**(1./n) * p.max_cryst_rate
    #
    phi = np.zeros(te.shape, dtype = te.dtype) # (rows, cols)
    phi[I] = 1 - np.exp(-(C_vent + lam*t_active)**n)
    return phi

def yield_stress(phi):
    # Mueller et al. (2010); Chevrel et al. (2013); FLOWGO
    #tau_star = 0.087 # depends on aspect ratio, technically: value of yield stress at phi = 0.293*phi_max
    #return tau_star * (np.maximum(1 - phi/p.phi_max, 1e-3)**-2 - 1)
    return 6500 * phi**2.85
    #return ne.evaluate('6500 * phi**2.85')

def slope_SAN(S, I):
    # Calculate abs(grad(S)): Steepest Adjacent Neighbor (SAN)
    # S.shape =  (rows, cols)
    #
    I_ij = I[1:-1,1:-1]
    S_ij = S[1:-1,1:-1][I_ij] # (n_active,)
    S_E = S[1:-1,2:][I_ij] # (n_active,)
    S_W = S[1:-1,0:-2][I_ij] # (n_active,)
    S_N = S[0:-2,1:-1][I_ij] # (n_active,)
    S_S = S[2:,1:-1][I_ij] # (n_active,)
    S_NE = S[0:-2,2:][I_ij] # (n_active,)
    S_SW = S[2:,0:-2][I_ij] # (n_active,)
    S_NW = S[0:-2,0:-2][I_ij] # (n_active,)
    S_SE = S[2:,2:][I_ij] # (n_active,)
    #
    dS_E = abs(S_E - S_ij) # (n_active,)
    dS_W = abs(S_ij - S_W) # (n_active,)
    dS_N = abs(S_N - S_ij) # (n_active,)
    dS_S = abs(S_ij - S_S) # (n_active,)
    dS_NE = abs(S_NE - S_ij) # (n_active,)
    dS_SW = abs(S_ij - S_SW) # (n_active,)
    dS_NW = abs(S_NW - S_ij) # (n_active,)
    dS_SE = abs(S_ij - S_SE) # (n_active,)
    #
    #diff_NE_SW = ne.evaluate('where(dS_NE > dS_SW, dS_NE, dS_SW)')
    #diff_NW_SE = ne.evaluate('where(dS_NW > dS_SE, dS_NW, dS_SE)')
    #slope_max_diag = ne.evaluate('where(diff_NE_SW > diff_NW_SE, diff_NE_SW, diff_NW_SE)') / np.sqrt(p.dx**2 + p.dy**2)
    #slope_E_W = ne.evaluate('where(dS_E > dS_W, dS_E, dS_W)') / p.dx
    #slope_N_S = ne.evaluate('where(dS_N > dS_S, dS_N, dS_S)') / p.dy
    #slope_max_card = ne.evaluate('where(slope_E_W > slope_N_S, slope_E_W, slope_N_S)')
    slope_max_diag = np.maximum.reduce([dS_NE, dS_SW, dS_NW, dS_SE]) / np.sqrt(p.dx**2 + p.dy**2) # (n_active,)
    slope_max_EW = np.maximum(dS_E, dS_W) / p.dx # (n_active,)
    slope_max_NS = np.maximum(dS_N, dS_S) / p.dy # (n_active,)
    #
    #return ne.evaluate('where(slope_max_card > slope_max_diag, slope_max_card, slope_max_diag)')
    return np.maximum.reduce([slope_max_EW, slope_max_NS, slope_max_diag])

def rheo_factor_bl_S(h, B, phi_S, cryst_core):
    #
    # --------------------------------------------------------------------------
    # Rheological Mobility factor
    # Contains all rheological information
    # related to gravity current diffusivity: K = R * h**3
    #
    Rs = np.zeros(h.shape, dtype = h.dtype) # (rows, cols)
    R = Rs.copy()
    fluidity_core = Rs.copy()
    abs_grad_S = Rs.copy()
    #
    valid = h > p.tiny_flow
    if np.any(valid):
        #
        grad_S_v = slope_SAN(B + h, valid) # (n_active,)
        phi_S_v = phi_S[valid] # (n_active,)
        h_v = h[valid] # (n_active,)
        cryst_core_v = cryst_core[valid] # (n_active,)
        #
        spec_weight = p.lava_density * p.g
        fluidity_core_v = fluidity(p.core_temperature, cryst_core_v)
        R_core = 0.3333 * spec_weight * fluidity_core_v
        #
        yield_strength = yield_stress(cryst_core_v) + phi_S_v*p.yield_strength_crust
        bingham_core = yield_strength / (spec_weight * h_v * grad_S_v + p.pos_eps) # Bingham Number w/ div-by-zero protection  # (n_active,)
        #
        term1 = np.maximum(1 - bingham_core, 0)**2
        term2 = np.maximum(phi_S_v - bingham_core, 0)**2
        f = term1 - term2
        fS_T_v = np.maximum(0, f)
        f_T_v = np.maximum(0, term1 - term2*phi_S_v + 0.5*bingham_core*f)
        #
        R[valid] = R_core * f_T_v
        Rs[valid] = 1.5 * R_core * fS_T_v
        fluidity_core[valid] = fluidity_core_v
        abs_grad_S[valid] = grad_S_v
    #
    return R, Rs, fluidity_core, abs_grad_S

def rheo_factor_bl_S_all_outputs(h, B, phi_S, cryst_core):
    #
    # --------------------------------------------------------------------------
    # Rheological Mobility factor
    # Contains all rheological information
    # related to gravity current diffusivity: K = R * h**3
    #
    Rs = np.zeros(h.shape, dtype = h.dtype) # (rows, cols)
    R = Rs.copy()
    fluidity_core = Rs.copy()
    abs_grad_S = Rs.copy()
    f_T = Rs.copy()
    #
    valid = h > p.tiny_flow
    if np.any(valid):
        #
        grad_S_v = slope_SAN(B + h, valid) # (n_active,)
        phi_S_v = phi_S[valid] # (n_active,)
        h_v = h[valid] # (n_active,)
        cryst_core_v = cryst_core[valid] # (n_active,)
        #
        spec_weight = p.lava_density * p.g
        fluidity_core_v = fluidity(p.core_temperature, cryst_core_v)
        R_core = 0.3333 * spec_weight * fluidity_core_v
        #
        yield_strength = yield_stress(cryst_core_v) + phi_S_v*p.yield_strength_crust
        bingham_core = yield_strength / (spec_weight * h_v * grad_S_v + p.pos_eps) # Bingham Number w/ div-by-zero protection  # (n_active,)
        #
        term1 = np.maximum(1 - bingham_core, 0)**2
        term2 = np.maximum(phi_S_v - bingham_core, 0)**2
        f = term1 - term2
        fS_T_v = np.maximum(0, f)
        f_T_v = np.maximum(0, term1 - term2*phi_S_v + 0.5*bingham_core*f)
        #
        R[valid] = R_core * f_T_v
        Rs[valid] = 1.5 * R_core * fS_T_v
        fluidity_core[valid] = fluidity_core_v
        abs_grad_S[valid] = grad_S_v
        f_T[valid] = f_T_v
    #
    return R, Rs, fluidity_core, abs_grad_S, f_T

def dBdt_constants(T_core):
    # Constants:
    tau_remelted = 3.35
    tau_transition = 0.64
    a_star = (T_core - p.ground_temperature) * d_log_visc_dT(T_core)
    #
    lam = therm.erfinv(1 - 1.38629 / a_star)
    temp1 = lam * 3.2 - 2.71
    temp2 = np.log(2.5 * lam)
    omega = 2./temp1 * (temp2 + np.sqrt(temp2**2 - temp1 / (lam * 1.8303)))
    A = -1./omega * np.log((np.exp(-omega * lam * 1.8303) - 1) / (1 - np.exp(0.8375*omega)))
    return a_star, lam, omega, A

def dBdt(t, ti, h, T_core, U_s, fluidity_core, q_n):
    # Constants:
    a_star, lam, omega, A = dBdt_constants(T_core)
    #
    out = np.zeros(h.shape)
    #
    valid = h > p.tiny_flow
    if np.any(valid):
        #
        h_v = h[valid]
        tb_v = t - ti[valid]
        fluidity_core_v = fluidity_core[valid]
        U_s_v = U_s[valid]
        #
        beta = 0.3197 / (a_star/0.6931 -  1) * (p.lava_density * p.lava_diffusivity * fluidity_core_v)**0.166667 #  == 0.3387 / (a_star/0.6931 -  1) * (1.4142 * kinematic_visc_core / p.lava_diffusivity)**-0.166667
        PeFo = U_s_v * tb_v / h_v
        X = beta**2 * PeFo
        E = np.exp(-omega * (2 * lam * X**0.5 - A + X))
        dBdt_cond = lam * (p.lava_diffusivity / tb_v)**0.5
        dBdt_conv = -beta * (U_s_v * p.lava_diffusivity / h_v)**0.5
        out[valid] = (dBdt_cond*E + dBdt_conv) / (E + 1)
        #
    #
    out[q_n > 0] = 0
    #
    return out






'''
notes on error for Forward Euler for dBdt:

max error occurs when dBdt ~ c * sqrt(1/t)

Error Anaylsis for:
dB/dt = sqrt(1/t)
B(dt/2) = 0
vs.
B_{n} = B_{n-1} + dt / sqrt(t_{n-1})
t_{0} = dt/2
B_{0} = 0

B(t_{n}) = sqrt(4*t_{n}) - sqrt(2*dt) = sqrt(4*(dt/2 + n*dt)) - sqrt(2*dt)
         = sqrt(dt) * sqrt(2) * (sqrt(2*n + 1) - 1)
B_{n} = sqrt(dt) * sum_{k=0}^{n-1} 1/sqrt(k + 1/2)

Global Error:
E_{n} = sqrt(dt) * ( sum_{k=0}^{n-1} 1/sqrt(k + 1/2) - sqrt(2)*(sqrt(2*n+1) -1) )

asympotics: n --> inf:
B_{n} ~ sqrt(dt) * (2*sqrt(n) + (sqrt(2)-1)*zeta(1/2)) + O(1/n**1.5)
E_{n} ~ sqrt(dt) * (sqrt(2) + (sqrt(2)-1)*zeta(1/2)) + O(1/sqrt(n))

lim n --> inf:
E_{n} --> sqrt(dt) * (sqrt(2) + (sqrt(2)-1)*zeta(1/2)) ~ 0.8093 * sqrt(dt)

'''



'''

def rheo_factor_bl(h, B, phi_B, phi_S, cryst_core, I):
    #
    # --------------------------------------------------------------------------
    # Rheological Mobility factor
    # Contains all rheological information
    # related to gravity current diffusivity: K = R * h**3
    #
    grad_S_active = slope_SAN(B + h, I) # (n_active,)
    #
    phi_B_active = phi_B[I] # (n_active,)
    phi_S_active = phi_S[I] # (n_active,)
    h_active = h[I] # (n_active,)
    cryst_core_active = cryst_core[I] # (n_active,)
    #
    spec_weight = p.lava_density * p.g
    R_core = 0.33333 * spec_weight / viscosity(p.core_temperature, cryst_core_active)

    bingham_core = yield_stress(cryst_core_active) / (spec_weight * h_active * grad_S_active + p.pos_eps) # Bingham Number w/ div-by-zero protection  # (n_active,)
    #
    term1 = np.maximum(1 - phi_B_active - bingham_core, 0)**2
    term2 = np.maximum(phi_S_active - bingham_core, 0)**2
    fS_T_active = np.maximum(0, term1 - term2)
    f_T_active = np.maximum(0, term1 * (1 - phi_B_active + 0.5*bingham_core) - term2 * (phi_S_active + 0.5*bingham_core))
    #f_T_active = np.maximum(0, ne.evaluate('term1 * (1 - phi_B_active + 0.5*bingham_core) - term2 * (phi_S_active + 0.5*bingham_core)'))
    #
    # Need to correct based on thickness threshold for valid cells in therm.thermal_BLs
    invalid_active = (h_active < p.tiny_flow)
    fS_T_active[invalid_active] = 0.
    f_T_active[invalid_active] = 0.
    #
    Rs = np.zeros(h.shape, dtype = h.dtype) # (rows, cols)
    R = Rs.copy()
    R[I] = R_core * f_T_active
    Rs[I] = 1.5 * R_core * fS_T_active
    #
    return R, Rs


'''

#

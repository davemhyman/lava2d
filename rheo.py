import numpy as np
from globals import params as p
import thermal as therm
################################################################################
# Rheological Model
################################################################################

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Specify Melt Viscosity Model
# must specify: viscosity(T), d_log_visc_dT(T), T_stiffened(T_core)

def melt_viscosity(T):
    # Dragoni (1989)
    a = d_log_visc_dT(p.core_temperature)
    return p.viscosity_vent * np.exp(a*(p.T_vent - T))

def d_log_visc_dT(T):
    # a_Dragoni = .04
    return .04

def T_stiffened(T_core):
    # assume all stiffening happens due to cooling of melt phase (no excess xtals)
    # ie: Temp at which:  visc_melt_stiff = 2*visc_melt_core
    return T_core - np.log(2) / d_log_visc_dT(T_core)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def viscosity(T, phi):
    return np.maximum(1 - phi/p.phi_max, 1e-3)**-2.5 * melt_viscosity(T)

def fluidity(T, phi):
    # reciprocal viscosity
    return np.maximum(1 - phi/p.phi_max, 0)**2.5 / melt_viscosity(T)

def cryst_avrami(t, te, I):
    # params, could be precalculated
    n = 4
    k = (p.max_cryst_rate/(n*p.phi_inf))**n * ((n-1.)/(n*np.exp(1)))**(1.-n)
    t_init = (np.log(p.phi_inf/(p.phi_inf - p.cryst_vent)) / k)**(1./n)
    #
    surface_age = t - te[I]
    t_bulk = 1.0 * surface_age # (n_active,) # core arrival in terms of surface transit time
    phi = np.zeros(te.shape, dtype = te.dtype) # (rows, cols)
    phi[I] = p.phi_inf * (1 - np.exp(-k*(t_init + t_bulk)**n))
    return phi

def yield_stress(phi):
    return 6500 * phi**2.85

def slope_SAN(S, I):
    # Calculate abs(grad(S)): Steepest Adjacent Neighbor (SAN)
    # S.shape =  (rows, cols)
    #
    I_ij = I[1:-1,1:-1]
    S_ij = S[1:-1,1:-1][I_ij] # (n_active,)
    #
    dS_E = abs(S_ij - S[1:-1,2:][I_ij]) # (n_active,)
    dS_W = abs(S_ij - S[1:-1,0:-2][I_ij]) # (n_active,)
    dS_N = abs(S_ij - S[0:-2,1:-1][I_ij]) # (n_active,)
    dS_S = abs(S_ij - S[2:,1:-1][I_ij]) # (n_active,)
    #
    slope_max_EW = np.maximum(dS_E, dS_W) / p.dx # (n_active,)
    slope_max_NS = np.maximum(dS_N, dS_S) / p.dy # (n_active,)
    slope_max = np.sqrt(slope_max_EW**2 + slope_max_NS**2)
    #
    return slope_max

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
    E_T = Rs.copy()
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
        T_stiff = T_stiffened(p.core_temperature) # temp at which visc = 2*visc_core
        eta_inf = therm.erfinv((T_stiff-p.atm_temperature)/(p.core_temperature-p.atm_temperature))
        phi_true_crust = phi_S_v * (0.88/eta_inf) # 0.88 from Hon et al., 1998
        yield_strength = yield_stress(cryst_core_v) + p.yield_strength_crust * phi_true_crust
        bingham_core = yield_strength / (spec_weight * h_v * grad_S_v + p.pos_eps) # Bingham Number w/ div-by-zero protection  # (n_active,)
        #
        term1 = np.maximum(1 - bingham_core, 0)**2
        term2 = np.maximum(phi_S_v - bingham_core, 0)**2
        f = term1 - term2
        ES_T_v = np.maximum(0, f)
        E_T_v = np.maximum(0, term1 - term2*phi_S_v + 0.5*bingham_core*f)
        #
        R[valid] = R_core * E_T_v
        Rs[valid] = 1.5 * R_core * ES_T_v
        fluidity_core[valid] = fluidity_core_v
        abs_grad_S[valid] = grad_S_v
        E_T[valid] = E_T_v
    #
    return R, Rs, fluidity_core, abs_grad_S, E_T

def dBdt_constants(T_core):
    # Constants:
    tau_remelted = 3.35
    tau_transition = 0.64
    a_star = (T_core - p.ground_temperature) * d_log_visc_dT(T_core)
    #
    lam = therm.erfinv(1 - 1.38629 / a_star)
    #
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
    valid = np.logical_and(h > p.tiny_flow, ti < t)
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
    lambda_substrate = 1.0 * lam
    out = lambda_substrate / lam * out
    #
    return out
#

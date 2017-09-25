import numpy as np
import scipy as sp
# check link below for kelvin functions
# http://docs.sympy.org/0.7.3/modules/mpmath/functions/bessel.html#ber
import mpmath as mp



# constansts
winding_factor = np.pi/(2*np.sqrt(3)) # for orthozyclic windings
mu_0                = 1.256 * 10**-6  # mag. feldkonstante [N/A^2]


# Funktionen:
# Radialgeschwindigkeit
def v_radial(omega, radius):
    #INPUT:  omega           -> angular speed [rad/s]
    #        radius          -> radius [m]
    #OUTPUT: v_radial        -> radial velocity [m/s]
    return omega * radius

#----------------------------------
#-----EIGENSCHAFTEN SPULE----------

# Innerer Widerstand Kupfer
def R_i(rho_mat, l_wire, A_wire):
    #INPUT:  rho_mat     -> sp. resistence of material [Ohm*m]
    #        l_wire      -> length of wire [m]
    #        A_wire      -> diameter of wire [m^2]
    #OUTPUT: R_i         -> inner resictance of wire [ohm]
    return rho_mat * l_wire / A_wire

# Resistance power due to Skin Effect (innerer Wirbelstrom)
def P_skin(I_peak,d_wire,f):
    #INPUT:  d      -> diameter of cable [m]
    #        f      -> frecuency of current [1/s]
    #OUTPUT: R_prox -> resictance due to prox-effect [ohm]
    # not sure if sigma of 1 is right
    sigma = 1
    delta = 1 / np.sqrt(np.pi * f * sigma * mu_0)
    xi = float(d_wire / (np.sqrt(2) * delta))
    R_DC = 4 / (sigma * np.pi * d_wire**2)

    kelvin =    (mp.ber(0, xi) * mp.bei(1, xi) - mp.ber(0, xi) * mp.ber(1, xi)) / (mp.ber(1, xi)**2 + mp.bei(1, xi)**2) \
              - (mp.bei(0, xi) * mp.ber(1, xi) + mp.bei(0, xi) * mp.bei(1, xi)) / (mp.ber(1, xi)**2 + mp.bei(1, xi)**2)

    F_R = xi / (4 * np.sqrt(2)) * kelvin

    return R_DC*F_R*I_peak


# Resistance power due to Proximity Effect (induzierter Wirbelstrom)
def P_prox(d_wire, f):
    #INPUT:  d_cable  -> diameter of cable [m]
    #        f        -> frecuency of current [1/s]
    #OUTPUT: R_prox   -> resictance due to prox-effect [ohm]
    sigma = 1
    delta = 1 / np.sqrt(np.pi * f * sigma * mu_0)
    xi = float(d_wire / (np.sqrt(2) * delta))
    R_DC = 4 / (sigma * np.pi * d_wire**2)

    kelvin =    (mp.ber(2, xi) * mp.ber(1, xi) + mp.ber(2, xi) * mp.bei(1, xi)) / (mp.ber(0, xi)**2 + mp.bei(0, xi)**2) \
              + (mp.bei(2, xi) * mp.bei(1, xi) + mp.bei(2, xi) * mp.ber(1, xi)) / (mp.ber(0, xi)**2 + mp.bei(0, xi)**2)

    G_R = xi * (np.pi * d_wire)**2 / (2 * np.sqrt(2)) * kelvin

    H_peak = 1 / (np.pi * d_wire)

    return R_DC * G_R * H_peak**2

def l_kabel(N, l_coil_turn ):
    #INPUT:  N              -> number of windings [1]
    #        l_coil_turn    -> length of wire per coil turn [m]
    #OUTPUT: l_kabel        -> lenght of cable [m]
    l_out = 1.5            #-> length outside of spool
    return (l_out +  N * l_coil_turn )

def A_kabel(d_kabel):
    #INPUT:  d_kabel        -> diameter of cable [m]
    #OUTPUT: A_kabel        -> cross-section area of cable [m^2]
    return ((d_kabel**2) / 4 * np.pi)

def x_mag(N, d_kabel, lagen):
    #INPUT:  N              -> number of windings [1]
    #        d_kabel        -> diameter of cable [m^2]
    #        lagen          -> num of layers [1]
    #OUTPUT: x_mag          -> average distance x to magnet [m]
    min_dist = 0.001        #-> min. distance is 1 mm
    return min_dist + (N / lagen * d_kabel) / 2

# B-Feld Berechnung
def B_peak(x_mag):
    #INPUT:  x_mag          -> average distance x to magnet [m]
    #OUTPUT: B_field        -> magnetic flow density [T]
    return 0.28 - 5 * x_mag

#---------------------------------------
#-------EIGENSCHAFTEN GENERATOR----------
# Effektivwert der Spannung
def e_eff(B_peak, l_eff, v_radial, N):
    #INPUT:  B_peak         -> peak value of magnetic flow density [T]
    #        l_eff          -> effective length of cable [m]
    #        v_radial       -> radial velocity [m/s]
    #        N              -> number of windings [1]
    #OUTPUT: e_eff          -> effective value of voltage [V]
    return B_peak* l_eff * v_radial * N * winding_factor/np.sqrt(2)

# induzierte Spannung
def e_indu(B, l_eff, v, N):
    #INPUT:  B              -> magnetic flow density [T]
    #        l_eff          -> effective length of cable [m]
    #        v_radial       -> radial velocity [m/s]
    #        N              -> number of windings [1]
    #OUTPUT: e_indu         -> induced voltage [V]
    return B * l_eff * v * N * winding_factor

# Induktivität
def L(N, r_outer, r_inner, l_eff, p):
    # FLäche = Außenkreisfläche - Innenkreisfläche * 60/360
    A = np.pi * (r_outer**2 - r_inner**2 ) * 60/360
    return mu_0 * N**2 * A / l_eff

# induzierter Strom
def i_indu(e, R_i, R_L, L, f):
    #INPUT:  e_indu  -> induced voltage [V]
    #        R_i     -> inner resistance of cable [ohm]
    #        R_L     -> load resistance [ohm]
    #        L       -> inductance of spool [H] 'Henry'
    #        f       -> frecuency of current [1/s]
    #OUTPUT: i_indu  -> induced amperage [A]
    return e / np.sqrt((R_i + R_L)**2 + (L * 2 * np.pi * f)**2)

# induzierter Strom
def i_eff(e_eff, R_i, R_L, L, f):
    #INPUT:  e_eff   -> effective voltage [V]
    #        R_i     -> inner resistance of cable [ohm]
    #        R_L     -> load resistance [ohm]
    #        L       -> inductance of spool [H] 'Henry'
    #        f       -> frecuency of current [1/s]
    #OUTPUT: i_indu  -> induced amperage [A]
    return e_eff / np.sqrt((R_i + R_L)**2 + (L * 2 * np.pi * f)**2)

# Frequenz
def f(n, p):
    return n * p

# Rotorleistung
def P_T(n, M_T):
    return n * M_T

# Lastleistung
def P_L(R_L, i):
    return R_L * i**2


def eta_wea(P_L, P_w):
    return P_L / P_w


def eta_g(P_L, P_T):
    return P_L / P_T

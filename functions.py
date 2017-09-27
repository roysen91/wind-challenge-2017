import numpy as np
import scipy as sp
# check link below for kelvin functions
# http://docs.sympy.org/0.7.3/modules/mpmath/functions/bessel.html#ber
import mpmath as mp



# constansts
winding_factor = np.pi/(2*np.sqrt(3)) # for orthozyclic windings
mu_0           = 1.256 * 10**-6  # mag. feldkonstante [N/A^2]
sigma_cu       = 58.5 # conductivity of copper [1/(Ohm*m)]


# Funktionen:
# Radialgeschwindigkeit
def v_radial(n, radius):
    #INPUT:  n           -> spool rotational speed [1/s]
    #        radius      -> radius [m]
    #OUTPUT: v_radial    -> radial velocity [m/s]
    return 2*np.pi*n * radius

#----------------------------------
#-----EIGENSCHAFTEN SPULE----------

# Innerer Widerstand 
def R_i(rho_mat, l_wire, A_wire):
    #INPUT:  rho_mat     -> sp. resistence of material [Ohm*m]
    #        l_wire      -> length of wire [m]
    #        A_wire      -> diameter of wire [mm^2]
    #OUTPUT: R_i         -> inner resictance of wire [ohm]
    return rho_mat * l_wire / A_wire

# Resistance power due to Skin Effect (innerer Wirbelstrom)
def P_skin(I_peak,d_wire,f,sigma=sigma_cu):
    #INPUT:  d      -> diameter of cable [mm]
    #        f      -> frecuency of current [1/s]
    #OUTPUT: R_prox -> resictance due to prox-effect [ohm]
    delta = 1 / np.sqrt(np.pi * f * sigma * mu_0)  # skin depth [m]
    xi = float(d_wire / (np.sqrt(2) * delta))
    R_DC = 4 / (sigma * np.pi * d_wire**2)

    kelvin =    (mp.ber(0, xi) * mp.bei(1, xi) - mp.ber(0, xi) * mp.ber(1, xi)) / (mp.ber(1, xi)**2 + mp.bei(1, xi)**2) \
              - (mp.bei(0, xi) * mp.ber(1, xi) + mp.bei(0, xi) * mp.bei(1, xi)) / (mp.ber(1, xi)**2 + mp.bei(1, xi)**2)

    F_R = xi / (4 * np.sqrt(2)) * kelvin

    return R_DC*F_R*I_peak**2


# Resistance power due to Proximity Effect (induzierter Wirbelstrom)
def P_prox(I_peak,d_wire, f, sigma = sigma_cu):
    #INPUT:  d_cable  -> diameter of cable [mm]
    #        f        -> frecuency of current [1/s]
    #OUTPUT: R_prox   -> resictance due to prox-effect [ohm]
    delta = 1 / np.sqrt(np.pi * f * sigma * mu_0)
    xi = float(d_wire / (np.sqrt(2) * delta))
    R_DC = 4 / (sigma * np.pi * d_wire**2)

    kelvin =    (mp.ber(2, xi) * mp.ber(1, xi) + mp.ber(2, xi) * mp.bei(1, xi)) / (mp.ber(0, xi)**2 + mp.bei(0, xi)**2) \
              + (mp.bei(2, xi) * mp.bei(1, xi) + mp.bei(2, xi) * mp.ber(1, xi)) / (mp.ber(0, xi)**2 + mp.bei(0, xi)**2)

    G_R = xi * (np.pi * d_wire)**2 / (2 * np.sqrt(2)) * kelvin

    H_peak = I_peak / (np.pi * d_wire)

    return R_DC * G_R * H_peak**2

def l_kabel(N, l_coil_turn ):
    #INPUT:  N              -> number of windings [1]
    #        l_coil_turn    -> length of wire per coil turn [m]
    #OUTPUT: l_kabel        -> lenght of cable [m]
    l_out = 1.5*10**-3     #-> length outside of spool [m]
    return (l_out +  N * l_coil_turn )

def A_kabel(d_kabel):
    #INPUT:  d_kabel        -> diameter of cable [mm]
    #OUTPUT: A_kabel        -> cross-section area of cable [mm^2]
    return (d_kabel/2)**2  * np.pi

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

    # B_peak from FEMM 0.65T in 1mm distance from magnet surface
    # B = 0.25 - 0.003*x_mag 
    return  0.25 - 0.003*x_mag 

#---------------------------------------
#-------EIGENSCHAFTEN GENERATOR----------

# Spannung
def Voltage(B, l_coil_eff, v, N):
    #INPUT:  B              -> magnetic flow density [T]
    #        l_coil_eff     -> effective length of wire per coil turn [m]
    #        v_radial       -> radial velocity [m/s]
    #        N              -> number of windings [1]
    #OUTPUT: e_indu         -> induced voltage [V]
    return B * l_coil_eff * v * N * winding_factor

# induzierter Strom
def Current(e, R_i, R_L, L, f):
    #INPUT:  e       -> voltage [V]
    #        R_i     -> inner resistance of cable [ohm]
    #        R_L     -> load resistance [ohm]
    #        L       -> inductance of spool [H] 'Henry'
    #        f       -> frecuency of current [1/s]
    #OUTPUT: i_peak  -> induced amperage [A]
    return e / np.sqrt((R_i + R_L)**2 + (L * 2 * np.pi * f)**2)

# Induktivität
def Inductance(N, r_outer, r_inner, l_eff, p):
    #INPUT:  N              -> number of windings [1]
    #        r_outer        -> outer radius of coil [m]
    #        r_inner        -> inner radius of coil [m]
    #        l_eff          ->  effective coil length per turn [m]
    #        p              -> number of pole pairs [1]
    #OUTPUT: L              -> inductivity [H] Henry
    # FLäche = Außenkreisfläche - Innenkreisfläche * 60/360
    A = np.pi * (r_outer**2 - r_inner**2 ) * 60/360
    return mu_0 * N**2 * A / 2*l_eff

# Frequenz
def f(n, p):
    #INPUT:  n    -> rotational speed of spool [1/s]
    #        p    -> number of pole pairs [1]
    #OUTPUT: f    -> frecuency of current [1/s]
    return n * p

# Rotorleistung
def P_T(n, M_T):
    return n * M_T

# Lastleistung
def P_L(R_L, i):
    return i**2 *R_L


def eta_wea(P_L, P_w):
    return P_L / P_w


def eta_g(P_L, P_T):
    return P_L / P_T

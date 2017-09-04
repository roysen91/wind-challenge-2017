import numpy as np
import scipy as sp


# Funktionen:

# Radialgeschwindigkeit
def v_radial(n, r_mag):
    return n * r_mag


#----------------------------------
#-----EIGENSCHAFTEN SPULE----------

# Innerer Widerstand Kupfer


def R_i(rho_cu, l_kabel, A_kabel):
    return rho_cu * l_kabel / A_kabel

# Proximity Effekt


def R_prox(d, f, mu_0):  # d = d_kabel
    sigma = 1
    delta = 1 / sqrt(np.pi * f * sigma * mu_0)
    xi = d / (sqrt(2) * delta)
    Rdc = 4 / (sigma * np.pi * d**2)
    kelvin = (sp.jv(2, xi) * sp.jv(1, xi) + sp.jv(2, xi) * sp.jn(1, xi)) / (sp.jv(0, xi)**2 + sp.jn(0, xi)**2) + (sp.jn(2, xi) * sp.jn(1, xi) + sp.jn(2, xi) * sp.jv(1, xi)) / (sp.jv(0, xi)**2 + sp.jn(0, xi)**2)
    Gr = xi * (np.pi * d)**2 / (2 * sqrt(2)) * kelvin
    H_amplitude = 1 / (np.pi * d_kabel)
    return Rdc * Gr * H_amplitude**2


def l_kabel(N, p, l_1, l_2, l_3):
    # 1.5m (Kabel außerhalb d. Spule), 2* weil PolPAAR, l_1 bis l_3 entspricht einer Wicklung
    return (1.5 + 2 * N * p * (l_1 + l_2 + l_3))

# Durchmesser d. Kabels


def A_kabel(d_kabel):
    return ((d_kabel**2) / 4 * np.pi)

# durchschnittlicher Abstand x zu Magnet


def x_mag(N, d_kabel, lagen):
    # der mittlere Abstand ist Minimalabstand v. 1 mm + Zahl der Wicklungen über Zahl der Lagen * des Durchmessers des Kabels und davon die Hälfte, da mittlerer Abstand
    return 0.001 + (N / lagen * d_kabel) / 2

# B-Feld Berechnung


def B_(x):
    return 0.28 - 5 * x

#---------------------------------------
#-------EIGENSCHAFTEN GENERATOR----------
# induzierte Spannung


def e_indu(B, l_eff, v, N, wickel):
    # Wicklungsgrad wickel entspricht pi/2*sqrt(3) = 0,9
    # bei orthozyklischer Wicklung
    return B * l_eff * v * N * wickel

# Induktivität


def L(N, l_eff, mu_0, p):
    # FLäche = Außenkreisfläche - Innenkreisfläche * 60/360
    A = (0.9**2 * np.pi - 0.45**2 * np.pi) * 120 / 360 / p
    return mu_0 * N**2 * A / l_eff

# induzierter Strom
def i_indu(e, R_V, R_L, L, f):
    return e / np.sqrt((R_V + R_L)**2 + (L * 2 * np.pi * f)**2)

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

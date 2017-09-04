import sys
sys.path.append(r'F:\\Programme\Anaconda\Lib')

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp


# Funktionen:

# Radialgeschwindigkeit
def _v_radial(n, r_mag):
    return n * r_mag


#----------------------------------
#-----EIGENSCHAFTEN SPULE----------

# Innerer Widerstand Kupfer


def _R_i(rho_cu, l_kabel, A_kabel):
    return rho_cu * l_kabel / A_kabel

# Proximity Effekt


def _R_prox(d, f, mu_0):  # d = d_kabel
    sigma = 1
    delta = 1 / sqrt(np.pi * f * sigma * mu_0)
    xi = d / (sqrt(2) * delta)
    Rdc = 4 / (sigma * np.pi * d**2)
    kelvin = (sp.jv(2, xi) * sp.jv(1, xi) + sp.jv(2, xi) * sp.jn(1, xi)) / (sp.jv(0, xi)**2 + sp.jn(0, xi)**2) + (sp.jn(2, xi) * sp.jn(1, xi) + sp.jn(2, xi) * sp.jv(1, xi)) / (sp.jv(0, xi)**2 + sp.jn(0, xi)**2)
    Gr = xi * (np.pi * d)**2 / (2 * sqrt(2)) * kelvin
    H_amplitude = 1 / (np.pi * d_kabel)
    return Rdc * Gr * H_amplitude**2


def _l_kabel(N, p, l_1, l_2, l_3):
    # 1.5m (Kabel außerhalb d. Spule), 2* weil PolPAAR, l_1 bis l_3 entspricht einer Wicklung
    return (1.5 + 2 * N * p * (l_1 + l_2 + l_3))

# Durchmesser d. Kabels


def _A_kabel(d_kabel):
    return ((d_kabel**2) / 4 * np.pi)

# durchschnittlicher Abstand x zu Magnet


def _x_mag(N, d_kabel, lagen):
    # der mittlere Abstand ist Minimalabstand v. 1 mm + Zahl der Wicklungen über Zahl der Lagen * des Durchmessers des Kabels und davon die Hälfte, da mittlerer Abstand
    return 0.001 + (N / lagen * d_kabel) / 2

# B-Feld Berechnung


def _B_(x):
    return 0.28 - 5 * x

#---------------------------------------
#-------EIGENSCHAFTEN GENERATOR----------
# induzierte Spannung


def _e_indu(B, l_eff, v, N, wickel):
    # Wicklungsgrad wickel entspricht pi/2*sqrt(3) = 0,9
    # bei orthozyklischer Wicklung
    return B * l_eff * v * N * wickel

# Induktivität


def _L(N, l_eff, mu_0, p):
    # FLäche = Außenkreisfläche - Innenkreisfläche * 60/360
    A = (0.9**2 * np.pi - 0.45**2 * np.pi) * 120 / 360 / p
    return mu_0 * N**2 * A / l_eff

# induzierter Strom


def _i_indu(e, R_V, R_L, L, f):
    return e / np.sqrt((R_V + R_L)**2 + (L * 2 * np.pi * f)**2)

# Frequenz


def _f(n, p):
    return n * p


# Rotorleistung


def _P_T(n, M_T):
    return n * M_T

# Lastleistung


def _P_L(R_L, i):
    return R_L * i**2


def _eta_wea(P_L, P_w):
    return P_L / P_w


def _eta_g(P_L, P_T):
    return P_L / P_T


# Stellgrößen (Werden nur verwendet, wenn diese konstant gehalten werden sollen):
p = 2  # Polpaarzahl
N = 20  # Wicklungen
Lagen = 10  # Zahl der Kabelschichten
d_kabel = 0.00175  # Kabeldurchmesser

# Variablen
n = 1485 / 60  # Umdrehungen d. Welle [1/s]
rho_cu = 17 * 10**-9  # sp. Widerstand Kupfer [Ohm+m]
mu_0 = 1.256 * 10**-6  # mag. feldkonstante
M_T = 6.924  # Moment Welle [Nm]
r_mag = 67.5 * 10**-3  # Radius Wellenmitte bis Magnetmitte
l_eff = 2 * 0.045

l_coil_aus = 0.9 * 2 * np.pi * 120 / 360 / p  # Außendurchmesser Spule
l_coil_inn = 0.45 * 2 * np.pi * 120 / 360 / p  # Innendurchmesser der Spule
P_w = 1.2 / 2 * 0.45**2 * np.pi * 10**3  # dichte luft = 1.2, rotorradius = 0.45, windgeschwindigkeit = 10
R_L = 10  # for schleife
wickel = 0.9  # Wicklungsdichte

x1 = _x_mag(N, d_kabel, Lagen)
B1 = _B_(x1)
e1 = _e_indu(B1, l_eff, _v_radial(n, r_mag), N, wickel)
l_kabel1 = _l_kabel(N, p, l_coil_inn, l_coil_aus, l_eff)
A1 = _A_kabel(d_kabel)
R_i1 = _R_i(rho_cu, l_kabel1, A1)
L1 = _L(N, l_eff, mu_0, p)
f1 = _f(n, p)
i1 = _i_indu(e1, R_i1, R_L, L1, f1)
P_v1 = R_i1 * i1**2
P_L1 = _P_L(R_L, i1)

print('Umdrehungen n =', n)
print('Dichte rho_cu =', rho_cu, 'kg/m^3')
print('Wellenmoment M_T =', M_T, 'Nm')
print('Radius Magnet r_mag =', r_mag * 1000, 'mm')
print('Effektive Länge l_eff =', l_eff * 1000, 'mm')
print('Länge Kabel l_kabel =', _l_kabel(N, p, l_coil_aus, l_coil_inn, l_eff), 'm')
print('Kabelfläche A =', round(A1 * 10**6, 5), 'mm^2')
print('Polpaarzahl =', p)
print('Mittelabstand Spule =', x1, 'm')
print('Magn. Flussdichte B =', B1, 'T')
print('Radialgeschwind. v =', _v_radial(n, r_mag), 'm/s')
print('induz. Spannung e =', e1, 'V')
print('induz. Strom i =', _i_indu(e1, R_i1, R_L, L1, f1), 'A')
print('Induktivität L =', L1, 'F')
print('Kupferwiderstand R_v =', R_i1, 'Ohm')
print('Frequenz f = ', f1, '1/s')
print('Windleistung P_L =', P_w, 'W')
print('Rotorleistung P_T =', _P_T(n, M_T), 'W')
print('Lastleistung P_L =', _P_L(R_L, i1), 'W')
print('Verlustleistung P_v =', P_v1, 'W')

# x_ar = np.linspace(0, 0.05, 51)
# plt.plot(x_ar, _B_(x_ar))
# plt.xlabel('Abstand [mm]')
# plt.ylabel('B [T]')
# plt.show()

N_ar = np.linspace(1, 300, 50)
plt.plot(N_ar, _R_i(rho_cu, _l_kabel(N_ar, p, l_coil_aus, l_coil_inn, l_eff), _A_kabel(d_kabel)), label='R_i')
plt.plot(N_ar, _l_kabel(N_ar, p, l_coil_aus, l_coil_inn, l_eff) / 100, label='Länge Kabel / 100')
plt.legend()
R_Last = [1, 5, 10, 20, 30, 40, 50]
# for i in range(0, 7):
#    it_legende = 'R_Last = ' + str(R_Last[i]) + ' Ohm'
#    it_e = _e_indu(_B_(_x_mag(N_ar, d_kabel, Lagen)), l_eff, _v_radial(n, r_mag), N_ar, wickel)
#    it_kabel = _l_kabel(N_ar, p, l_coil_aus, l_coil_inn, l_eff)
#    it_Ri = _R_i(rho_cu, it_kabel, _A_kabel(d_kabel))
#    it_i = _i_indu(it_e, it_Ri, R_Last[i], _L(N_ar, l_eff, mu_0, p), _f(n, p))
#    plt.plot(N_ar, _eta_g(_P_L(R_Last[i], it_i), P_w), label=it_legende)
#    plt.legend()
d_kabel_ar = np.linspace(0.001, 0.0025, 7)
# for i in range(0, 7):
#    it_legende = 'd_kabel = ' + str(round(d_kabel_ar[i] * 1000, 2)) + ' mm'
#    it_e = _e_indu(_B_(_x_mag(N_ar, d_kabel_ar[i], Lagen)), l_eff, _v_radial(n, r_mag), N_ar, wickel)
#    it_kabel = _l_kabel(N_ar, p, l_coil_aus, l_coil_inn, l_eff)
#    it_Ri = _R_i(rho_cu, it_kabel, _A_kabel(d_kabel_ar[i]))
#    it_i = _i_indu(it_e, it_Ri, R_Last[6], _L(N_ar, l_eff, mu_0, p), _f(n, p))
#    plt.plot(N_ar, _eta_g(_P_L(R_Last[6], it_i), P_w), label=it_legende)
#    plt.legend()
plt.xlabel('Wicklungen')
plt.ylabel('Ohm / m')
plt.show()


import sys
sys.path.append(r'F:\\Programme\Anaconda\Lib')

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import functions as fcn


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

x1 = fcn.x_mag(N, d_kabel, Lagen)
B1 = fcn.B_field(x1)
e1 = fcn.e_indu(B1, l_eff, fcn.v_radial(n, r_mag), N)
l_kabel1 = fcn.l_kabel(N, p, l_coil_inn, l_coil_aus, l_eff)
A1 = fcn.A_kabel(d_kabel)
R_i1 = fcn.R_i(rho_cu, l_kabel1, A1)
L1 = fcn.L(N, l_eff, p)
f1 = fcn.f(n, p)
i1 = fcn.i_indu(e1, R_i1, R_L, L1, f1)
P_v1 = R_i1 * i1**2
P_L1 = fcn.P_L(R_L, i1)

print('Umdrehungen n =', n)
print('Dichte rho_cu =', rho_cu, 'kg/m^3')
print('Wellenmoment M_T =', M_T, 'Nm')
print('Radius Magnet r_mag =', r_mag * 1000, 'mm')
print('Effektive Länge l_eff =', l_eff * 1000, 'mm')
print('Länge Kabel l_kabel =', fcn.l_kabel(N, p, l_coil_aus, l_coil_inn, l_eff), 'm')
print('Kabelfläche A =', round(A1 * 10**6, 5), 'mm^2')
print('Polpaarzahl =', p)
print('Mittelabstand Spule =', x1, 'm')
print('Magn. Flussdichte B =', B1, 'T')
print('Radialgeschwind. v =', fcn.v_radial(n, r_mag), 'm/s')
print('induz. Spannung e =', e1, 'V')
print('induz. Strom i =', fcn.i_indu(e1, R_i1, R_L, L1, f1), 'A')
print('Induktivität L =', L1, 'F')
print('Kupferwiderstand R_v =', R_i1, 'Ohm')
print('Frequenz f = ', f1, '1/s')
print('Windleistung P_L =', P_w, 'W')
print('Rotorleistung P_T =', fcn.P_T(n, M_T), 'W')
print('Lastleistung P_L =', fcn.P_L(R_L, i1), 'W')
print('Verlustleistung P_v =', P_v1, 'W')

# x_ar = np.linspace(0, 0.05, 51)
# plt.plot(x_ar, _B_field(x_ar))
# plt.xlabel('Abstand [mm]')
# plt.ylabel('B [T]')
# plt.show()

N_ar = np.linspace(1, 300, 50)
plt.plot(N_ar, fcn.R_i(rho_cu, fcn.l_kabel(N_ar, p, l_coil_aus, l_coil_inn, l_eff), fcn.A_kabel(d_kabel)), label='R_i')
plt.plot(N_ar, fcn.l_kabel(N_ar, p, l_coil_aus, l_coil_inn, l_eff) / 100, label='Länge Kabel / 100')
plt.legend()
R_Last = [1, 5, 10, 20, 30, 40, 50]
# for i in range(0, 7):
#    it_legende = 'R_Last = ' + str(R_Last[i]) + ' Ohm'
#    it_e = _e_indu(_B_field(_x_mag(N_ar, d_kabel, Lagen)), l_eff, _v_radial(n, r_mag), N_ar)
#    it_kabel = _l_kabel(N_ar, p, l_coil_aus, l_coil_inn, l_eff)
#    it_Ri = _R_i(rho_cu, it_kabel, _A_kabel(d_kabel))
#    it_i = _i_indu(it_e, it_Ri, R_Last[i], _L(N_ar, l_eff, mu_0, p), _f(n, p))
#    plt.plot(N_ar, _eta_g(_P_L(R_Last[i], it_i), P_w), label=it_legende)
#    plt.legend()
d_kabel_ar = np.linspace(0.001, 0.0025, 7)
# for i in range(0, 7):
#    it_legende = 'd_kabel = ' + str(round(d_kabel_ar[i] * 1000, 2)) + ' mm'
#    it_e = _e_indu(_B_field(_x_mag(N_ar, d_kabel_ar[i], Lagen)), l_eff, _v_radial(n, r_mag), N_ar)
#    it_kabel = _l_kabel(N_ar, p, l_coil_aus, l_coil_inn, l_eff)
#    it_Ri = _R_i(rho_cu, it_kabel, _A_kabel(d_kabel_ar[i]))
#    it_i = _i_indu(it_e, it_Ri, R_Last[6], _L(N_ar, l_eff, mu_0, p), _f(n, p))
#    plt.plot(N_ar, _eta_g(_P_L(R_Last[6], it_i), P_w), label=it_legende)
#    plt.legend()
plt.xlabel('Wicklungen')
plt.ylabel('Ohm / m')
plt.show()


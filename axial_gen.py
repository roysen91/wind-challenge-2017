import numpy as np

# Constants
num_pole_pairs 		= 2  # Polpaarzahl
num_coils 			= 4 # number of coils used [1]

# Variable generator parameter
N 					= 7  # Wicklungen
lagen 				= 4  # Zahl der Kabelschichten
d_kabel		 		= np.array([1.2, 1.5 ]) #[mm]
rho_cu 				= np.array([15.11 * 10**-3, 9.756 * 10**-3])  # sp. Widerstand Kupfer [Ohm/m] see Generator v10.pdf

# Radialgenerator
r_inner 			= 45.5*10**-3 # inner radius of magnets [mm]
r_outer 			= 90.5*10**-3 # outer radius of magnets [mm]
r_mag 				= (r_inner+r_outer)/2  # Radius Wellenmitte bis Magnetmitte [m]

l_coil_eff 			= 2 * 0.045		# effective coil length per turn [m]
l_coil_outer 		= (r_outer * 2 * np.pi / 360) * 60  # outer coil length per turn [m]
l_coil_inner 		= (r_inner * 2 * np.pi / 360) * 60# inner coil length per turn [m]
l_coil_turn 		= l_coil_eff + l_coil_outer + l_coil_inner

# Turbine data
M_T_N 	= 1.11 	# nominal turbine torque [Nm]
rps_T_N = 24.76 # nominal turbine rpm [1/s]


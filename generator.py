import numpy as np
import scipy as sp
# check link below for kelvin functions
# http://docs.sympy.org/0.7.3/modules/mpmath/functions/bessel.html#ber
import mpmath as mp

mu_0           = 1.256 * 10**-6  # mag. feldkonstante [N/A^2]


class Generator:
	def __init__(self,design = 'axial',pole_pairs=2,coils=4):
		self.design 		= design
		self.num_pole_pairs = pole_pairs
		self.num_coils 		= coils
		self.rot_speed 		= 24.76 # nominal turbine rpm [1/s]
		self.R_L			= 20 

		if self.design == 'radial':
			# Rotor geometry
			self.rotor_r_inner 		= 35*10**-3 # inner radius of magnets [mm]
			self.rotor_r_outer 		= 45*10**-3 # outer radius of magnets [mm]
			self.r_magnet 			= (self.rotor_r_inner+self.rotor_r_outer)/2  # radius to mid magnet [m]
			# stator geometry
			self.stator_r_inner 	= 47*10**-3 # inner radius [mm]
			self.stator_r_outer 	= 50*10**-3 # outer radius [mm]
			self.dist_rot_stat 	= 0.001+self.stator_r_outer-self.stator_r_inner  #-> min. distance is 1 mm + stator material thickness
			# coil geometry
			l_coil_eff 			= 2 * 0.120		# effective coil length per turn [m]
			l_coil_outer 		= (self.stator_r_outer * 2 * np.pi / 360) * 60  # outer coil length per turn [m]
			l_coil_inner 		= l_coil_outer  # inner coil length per turn [m]
			l_coil_turn 		= l_coil_eff + l_coil_outer + l_coil_inner
			area_coil 			= l_coil_eff * 2*np.pi*self.stator_r_outer *(60/360)

		else:
			# Rotor geometry
			self.rotor_r_inner 		= 45.5*10**-3 # inner radius of magnets [mm]
			self.rotor_r_outer 		= 90.5*10**-3 # outer radius of magnets [mm]
			self.r_magnet 			= (self.rotor_r_inner+self.rotor_r_outer)/2  # radius to mid magnet [m]
			# stator geometry
			self.stator_r_inner 	= self.rotor_r_inner # inner radius [mm]
			self.stator_r_outer 	= self.rotor_r_outer # outer radius [mm]
			
			self.dist_rot_stat 	= 0.001        #-> min. distance is 1 mm
			# coil geometry
			l_coil_eff 			= 2 * 0.045		# effective coil length per turn [m]
			l_coil_outer 		= (self.stator_r_outer * 2 * np.pi / 360) * 60  # outer coil length per turn [m]
			l_coil_inner 		= (self.stator_r_inner * 2 * np.pi / 360) * 60# inner coil length per turn [m]
			l_coil_turn 		= l_coil_eff + l_coil_outer + l_coil_inner
			area_coil 			= np.pi * (self.stator_r_outer**2 - self.stator_r_inner**2 ) * 60/360
			 
		self.coil = Coil(l_coil_eff,l_coil_turn,area_coil)

		self.compute()

	def compute(self):
		self.coil.compute()
		self.omega 		= self.rot_speed*2*np.pi
		self.f_el 		= self.frequency()
		self.x_mag 		= self.dist_mag()
		self.B			= self.mag_flow_dens()
		self.v_magnet 	= self.v_radial(self.r_magnet)
		self.V 			= self.voltage()
		self.I 			= self.current()

		self.P_L 		= self.R_L*self.I**2
		self.P_skin 	= self.power_skin()
		self.P_prox 	= self.power_prox()
		self.P_wire 	= self.num_coils*self.coil.R*self.I**2
		self.P_loss		= self.P_skin+self.P_prox+self.P_wire

		self.M 			= (self.P_L+self.P_loss)/self.omega

		self.eff 		= 1/(1+self.P_loss/(self.V*self.I))


	def v_radial(self,radius):
		#INPUT:  radius      -> radius [m]
	    #OUTPUT: v_radial    -> radial velocity [m/s]
	    return 2*np.pi*self.rot_speed * radius

	def frequency(self):
	    #OUTPUT: f    -> frecuency of current [1/s]
	    return self.rot_speed * self.num_pole_pairs
		
	def dist_mag(self):
	    #OUTPUT: x_mag          -> distance from magnet surface to mid of coil[m]
	    return self.dist_rot_stat + self.coil.layers * self.coil.d_wire

	def mag_flow_dens(self):
	    #OUTPUT: B_field        -> magnetic flow density [T]
	    # B_peak from FEMM 0.65T in 1mm distance from magnet surface
	    # B = 0.25 - 0.003*x_mag 
	    return  0.25 - 0.003*self.x_mag

	def voltage(self):
	    #OUTPUT: e_indu         -> induced voltage [V]
	    return self.B * self.coil.l_eff * self.v_magnet * self.coil.windings * self.coil.winding_factor

	def current(self):
	    #OUTPUT: i_peak  -> induced amperage [A]
	    return self.num_coils*self.V / np.sqrt((self.num_coils*self.coil.R + self.R_L)**2 + (self.num_coils*self.coil.L * 2 * np.pi * self.f_el)**2)

	def power_skin(self):
	    #OUTPUT: R_prox -> resictance due to prox-effect [ohm]
	    delta 	= 1 / np.sqrt(np.pi * self.f_el * self.coil.sigma_mat * mu_0)  # skin depth [m]
	    xi 		= float(self.coil.d_wire / (np.sqrt(2) * delta))
	    R_DC 	= 4 / (self.coil.sigma_mat * np.pi * self.coil.d_wire**2)

	    kelvin 	=   (mp.ber(0, xi) * mp.bei(1, xi) - mp.ber(0, xi) * mp.ber(1, xi)) / (mp.ber(1, xi)**2 + mp.bei(1, xi)**2) \
	              - (mp.bei(0, xi) * mp.ber(1, xi) + mp.bei(0, xi) * mp.bei(1, xi)) / (mp.ber(1, xi)**2 + mp.bei(1, xi)**2)

	    F_R 	= xi / (4 * np.sqrt(2)) * kelvin

	    return R_DC*F_R*self.I**2

	def power_prox(self):
	    #OUTPUT: R_prox   -> resictance due to prox-effect [ohm]
	    delta 	= 1 / np.sqrt(np.pi * self.f_el * self.coil.sigma_mat * mu_0)
	    xi 		= float(self.coil.d_wire / (np.sqrt(2) * delta))
	    R_DC 	= 4 / (self.coil.sigma_mat * np.pi * self.coil.d_wire**2)

	    kelvin 	=    (mp.ber(2, xi) * mp.ber(1, xi) + mp.ber(2, xi) * mp.bei(1, xi)) / (mp.ber(0, xi)**2 + mp.bei(0, xi)**2) \
	              + (mp.bei(2, xi) * mp.bei(1, xi) + mp.bei(2, xi) * mp.ber(1, xi)) / (mp.ber(0, xi)**2 + mp.bei(0, xi)**2)

	    G_R 	= xi * (np.pi * self.coil.d_wire)**2 / (2 * np.sqrt(2)) * kelvin

	    H_peak 	= self.I / (np.pi * self.coil.d_wire)

	    return R_DC * G_R * H_peak**2

class Coil:
	def __init__(self,l_eff,l_turn,area,l_out=1.5*10**-3,N=8,layers=5,d_wire=1.5,rho_mat=9.756 * 10**-3):
		self.l_eff 			= l_eff
		self.l_turn 		= l_turn
		self.area 			= area
		self.windings 		= N
		self.layers 		= layers
		self.d_wire			= d_wire
		self.l_out 			= l_out    #-> length outside of spool [m]
		self.winding_factor = np.pi/(2*np.sqrt(3)) # for orthozyclic windings
		self.rho_mat 		= rho_mat # sp. Widerstand Kupfer [Ohm/m] see Generator v10.pdf

		self.compute()

	def compute(self):
		self.sigma_mat      = 1/self.rho_mat # conductivity of copper [1/(Ohm*m)]
		self.l_wire 		= self.length_wire()
		self.A_wire 		= self.area_wire()
		self.R 				= self.resistance()
		self.L 				= self.inductance()

	def length_wire(self):
	    #OUTPUT: l_wire        -> lenght of wire [m]
	    return (self.l_out +  self.windings * self.l_turn )

	def area_wire(self):
	    #OUTPUT: area_wire      -> cross-section area of wire [mm^2]
	    return (self.d_wire/2)**2  * np.pi

	def resistance(self):
	    #OUTPUT: R_i         -> inner resictance of wire [ohm]
	    return self.rho_mat * self.l_wire / self.A_wire

	def inductance(self):
	    #OUTPUT: L              -> inductivity [H] Henry
	    # FLäche = Außenkreisfläche - Innenkreisfläche * 60/360
	    return mu_0 * self.windings**2 * self.area / 2*self.l_eff


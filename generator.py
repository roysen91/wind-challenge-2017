import numpy as np
import scipy as sp
# check link below for kelvin functions
# http://docs.sympy.org/0.7.3/modules/mpmath/functions/bessel.html#ber
import mpmath as mp

mu_0 			= 1.256 * 10**-6  # mag. feldkonstante [N/A^2]
rho_cu 			= 0.0171 # [Ohm*mm2/m] spezifischer Widerstand Kupfer


class Generator:
	def __init__(self,design = 'axial',pole_pairs=2,coils=4):
		self.design 		= design
		self.num_pole_pairs = pole_pairs
		self.num_coils 		= coils
		self.rot_speed 		= 24.76 # nominal turbine rpm [1/s]
		self.M_T 			= 1.1 # [Nm] form turbine 
		self.R_L			= 5 #[Ohm]

		if self.design == 'radial':
			self.b_avg 				= np.loadtxt('b_avg_radial_gen.dat') # load B-avg from file 
			self.angle_magnet 		= 70 # angle of magnets
			self.angle_space 		= 180/self.num_pole_pairs-self.angle_magnet # angle of space between magnets
			self.angle_coil 		= 20 # angle of coils
			self.angle_coil_space 	= (360/self.num_coils)-self.angle_coil # angle of space between coils
			# Rotor geometry
			self.rotor_r_inner 		= 35*10**-3 # inner radius of magnets [m]
			self.rotor_r_outer 		= 45*10**-3 # outer radius of magnets [m]
			self.r_magnet 			= (self.rotor_r_inner+self.rotor_r_outer)/2  # radius to mid magnet [m]
			# stator geometry
			self.stator_r_inner 	= 47*10**-3 # inner radius [m]
			self.stator_r_outer 	= 50*10**-3 # outer radius [m]
			self.dist_rot_stat 		= 0.002+self.stator_r_outer-self.stator_r_inner  #-> min. distance is 2 mm + stator material thickness
			# to end at 50% of magnet length where magnetic field is strongest
			l_coil_eff 				= 0.120		# effective coil length per side [m]
			l_coil_outer 			= (self.stator_r_outer * 2*np.pi/360)*(self.angle_coil+self.angle_coil_space)  # outer coil length per turn [m]
			l_coil_inner 			= (self.stator_r_outer * 2*np.pi/360)*(self.angle_coil+self.angle_coil_space)  # inner coil length per turn [m]
			l_coil_space 			= 2*self.stator_r_outer*np.pi*(self.angle_coil_space/360)
			# max width of layers above magnet surface
			max_coil_width = (2*np.pi*self.stator_r_outer*self.angle_magnet/360)*0.8

		else:
			self.b_avg 				= np.loadtxt('b_avg_axial_gen.dat') # load B-avg from file 
			self.angle_magnet 		= 60 # angle of magnets
			self.angle_space 		= 180/self.num_pole_pairs-self.angle_magnet # angle of space between magnets
			# Rotor geometry
			self.rotor_r_inner 		= 45.5*10**-3 # inner radius of magnets [m]
			self.rotor_r_outer 		= 90.5*10**-3 # outer radius of magnets [m]
			self.r_magnet 			= (self.rotor_r_inner+self.rotor_r_outer)/2  # radius to mid magnet [m]
			# stator geometry
			self.stator_r_inner 	= self.rotor_r_inner # inner radius [m]
			self.stator_r_outer 	= self.rotor_r_outer # outer radius [m]
			
			self.dist_rot_stat 		= 0.001        #-> min. distance is 1 mm
			# coil geometry
			l_coil_eff 				= 0.045		# effective coil length per side [m]
			# to end at 50% of magnet length where magnetic field is strongest
			l_coil_outer 			= (self.stator_r_outer * 2*np.pi/360)*(self.angle_magnet+self.angle_space) # outer coil length per turn [m]
			l_coil_inner 			= (self.stator_r_inner * 2*np.pi/360)*(self.angle_magnet+self.angle_space) # inner coil length per turn [m]
			l_coil_space 			= (self.rotor_r_inner+self.rotor_r_outer)*np.pi*(self.angle_space/360)
			# max width of layers above magnet surface
			max_coil_width = l_coil_inner*0.8


			 
		self.coil = Coil(l_coil_outer,l_coil_inner,l_coil_eff,l_coil_space,max_coil_width,self.num_coils)

		self.compute()

	def compute(self):
		self.coil.compute()
		self.l_wire 	= self.num_coils*self.coil.l_wire
		self.R_i 		= self.num_coils*self.coil.R
		self.omega 		= self.rot_speed*2*np.pi
		self.f_el 		= self.frequency()
		self.x_mag_l,self.x_mag_u 		= self.dist_mag()
		self.B_eff		= self.mag_flow_dens()
		self.H_peak 	= self.mag_field_strength()*np.sqrt(2)
		self.v_magnet 	= self.v_radial(self.r_magnet)
		self.V_eff 		= self.num_coils*self.voltage()
		self.wave_resist_coeff()
		

		self.P_in 		= self.omega*self.M_T# P_in is from Turbine
		self.I_eff 		= self.current()
		self.R_L 		= self.load_resist()

		self.P_skin 	= 2*self.F_R*self.num_coils*self.coil.R*self.I_eff**2
		self.P_prox 	= self.G_R*self.num_coils*self.coil.R*self.H_peak**2
		self.P_loss		= self.P_skin+self.P_prox
		self.P_L 		= self.R_L*self.I_eff**2

		if self.I_eff < 0 or self.R_L < 1 or self.R_L > 50:
			self.eff 	= 0
			self.M 		= 0
		else:
			self.eff 	= 1-self.P_loss/self.P_in
			self.M 		= (self.P_L+self.P_loss)/self.omega

	def v_radial(self,radius):
		#INPUT:  radius      -> radius [m]
	    #OUTPUT: v_radial    -> radial velocity [m/s]
	    return 2*np.pi*self.rot_speed * radius

	def frequency(self):
	    #OUTPUT: f    -> frecuency of current [1/s]
	    return self.rot_speed * self.num_pole_pairs
		
	def dist_mag(self):
	    #OUTPUT: x_mag          -> distance from magnet surface to mid of coil[m]
	    x_mag_l = self.dist_rot_stat + self.coil.d_wire/2
	    x_mag_u = self.dist_rot_stat + self.coil.d_wire*self.coil.layer_y-self.coil.d_wire/2
	    return x_mag_l,x_mag_u

	def mag_flow_dens(self):
	    #OUTPUT: B_eff        -> magnetic flow density [T]
	    B = (abs(np.interp(self.x_mag_l,self.b_avg[:,0],self.b_avg[:,1])+abs(np.interp(self.x_mag_u,self.b_avg[:,0],self.b_avg[:,1]))))/2
	    return B/np.sqrt(2)

	def voltage(self):
	    #OUTPUT: e_indu         -> induced voltage [V]
	    return self.B_eff * self.coil.l_eff * self.v_magnet * self.coil.windings * self.coil.winding_factor*0.9

	def current(self):
	    #OUTPUT: i_eff  -> induced amperage [A]
	    #return self.V_eff / np.sqrt((self.num_coils*self.coil.R + self.R_L)**2 + (self.num_coils*self.coil.L * 2 * np.pi * self.f_el)**2)
	    
	    p = self.V_eff/(2*self.R_i*self.F_R-self.R_i)
	    q = +(self.R_i*self.G_R*self.H_peak**2-self.P_in)/(2*self.R_i*self.F_R-self.R_i)

	    self.I_1 = (-p/2)+np.sqrt((p/2)**2-q)
	    self.I_2 = (-p/2)-np.sqrt((p/2)**2-q)
	    
	    return self.I_1

	def load_resist(self):
		#OUTPUT: R_load -> resictance due to load [ohm]
		return np.sqrt((self.V_eff/self.I_eff)**2-(self.num_coils*self.coil.L * 2 * np.pi * self.f_el)**2)-self.R_i

	def wave_resist_coeff(self):
	    #OUTPUT: F_R -> coeff for skin effekt [ohm]
		#		 G_R -> coeff for prox effekt [ohm]

		delta 	= 1 / np.sqrt(np.pi * self.f_el * (1/(self.coil.rho_mat*10**-6)) * mu_0)  # skin depth [m]
		xi 		= self.coil.d_wire / (np.sqrt(2) * delta)

		kelvin_skin 	= (mp.ber(0, xi) * mp.bei(1, xi) - mp.ber(0, xi) * mp.ber(1, xi)) / (mp.ber(1, xi)**2 + mp.bei(1, xi)**2) \
						- (mp.bei(0, xi) * mp.ber(1, xi) + mp.bei(0, xi) * mp.bei(1, xi)) / (mp.ber(1, xi)**2 + mp.bei(1, xi)**2)

		kelvin_prox 	= (mp.ber(2, xi) * mp.ber(1, xi) + mp.ber(2, xi) * mp.bei(1, xi)) / (mp.ber(0, xi)**2 + mp.bei(0, xi)**2) \
						+ (mp.bei(2, xi) * mp.bei(1, xi) - mp.bei(2, xi) * mp.ber(1, xi)) / (mp.ber(0, xi)**2 + mp.bei(0, xi)**2)

		self.F_R 		= xi / (4 * np.sqrt(2)) * kelvin_skin
		self.G_R 		= -(xi * (np.pi * self.coil.d_wire)**2 / (2 * np.sqrt(2))) * kelvin_prox

	def mag_field_strength(self):
		#OUTPUT: H -> magnetic field strength [A/m]
		return self.B_eff/mu_0

	def __str__(self):
		return 'l_wire 	= {} m\nx_layer = {}\ny_layer = {} \nR_i 	= {} Ohm\nB_eff 	= {} T\nV_eff 	= {} V\nI_eff 	= {} A\nR_L 	= {} Ohm\neta 	= {} %\n'\
				.format(np.round(self.l_wire,2),np.round(self.coil.layer_x,2),np.round(self.coil.layer_y,2),np.round(self.R_i,2),np.round(self.B_eff,2),np.round(self.V_eff,2),np.round(float(self.I_eff),2),np.round(float(self.R_L),2),np.round(float(self.eff)*100,2))



class Coil:
	def __init__(self,l_coil_outer,l_coil_inner,l_eff,l_space,max_coil_width,num_coils,l_out=1.5*10**-3,N=30,d_wire=0.001,rho_mat=rho_cu):
		self.max_coil_width = max_coil_width
		self.num_coils 		= num_coils
		self.min_coil_dist 	= 0.001 # [m] minimal distance between to coils
		self.l_eff			= 2*l_eff
		self.l_inner 		= l_coil_inner # [m] inner radius
		self.l_outer 		= l_coil_outer # [m] outer radius
		self.set_d_wire(d_wire) # [mm]
		self.set_windings(N)
		self.l_out 			= l_out    #-> length outside of spool [m]
		self.winding_factor = 1#((self.d_wire/2)**2*np.pi)/(1.32*10**-3)**2 # for rectengular winding with 0.3mm isolation thickness
		self.rho_mat 		= rho_mat # sp. Widerstand Kupfer [Ohm*mm^2/m] see Generator v10.pdf
		self.l_space 		= l_space # average distance between two coil cores

		self.compute()

	def compute(self):
		self.sigma_mat      = 1/self.rho_mat # conductivity of copper [1/(Ohm*m)]
		self.l_wire 		= self.length_wire() 	# [m]
		self.A_wire 		= self.area_wire()		# [mm^2]
		self.R 				= self.resistance()		# [Ohm]
		self.L 				= self.inductance()		# [H]

	def set_d_wire(self,d_wire,d_iso=0.07*10**-3):
		self.d_wire = d_wire # [mm]
		self.d_iso = d_iso # [mm]

	def set_windings(self,N):
		if self.num_coils<=2:
			self.layer_x_max = int(self.max_coil_width/(self.d_wire+self.d_iso)) # layers along magnet surface
		else:
			self.layer_x_max = int(self.max_coil_width/(2*(self.d_wire+self.d_iso))) # layers along magnet surface

		self.windings 		= N

		if self.windings<self.layer_x_max:
			self.layer_x 	= self.windings
			self.layer_y 	= 1
		else:
			self.layer_x 	= self.layer_x_max
			self.layer_y 	= int(self.windings/self.layer_x)+1
		self.l_turn 		= (self.l_outer+self.layer_x*(self.d_wire+self.d_iso))+(self.l_inner+self.layer_y*(self.d_wire+self.d_iso))+self.l_eff
		self.area 			= self.windings*self.d_wire#area 	 # [m^2]

	def check_geometry(self):
		# check if for number of windings coils have at least a minimum distance
		if self.num_coils>2:
			if (self.max_coil_width-2*self.layer_x*self.d_wire) <= self.min_coil_dist:
				return False
			else:
				return True
		else:
			if (self.layer_x*self.d_wire) <= self.max_coil_width:
				return False
			else:
				return True
				


	def length_wire(self):
	    #OUTPUT: l_wire        -> lenght of wire [m]
	    return (2*self.l_out +  self.windings * self.l_turn )*0.84

	def area_wire(self):
	    #OUTPUT: area_wire      -> cross-section area of wire [mm^2]
	    return ((self.d_wire*10**3)/2)**2  * np.pi

	def resistance(self):
	    #OUTPUT: R_i         -> inner resictance of wire [ohm]
	    return self.rho_mat * self.l_wire / self.A_wire

	def inductance(self):
	    #OUTPUT: L              -> inductivity [H] Henry
	    r_w = 45/2 #
	    #return mu_0 * self.windings**2 * self.area / self.l_eff
	    return mu_0 * self.windings**2 * self.area / (self.l_eff+2*r_w/2.2)


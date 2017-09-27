import numpy as np
import scipy as sp
# check link below for kelvin functions
# http://docs.sympy.org/0.7.3/modules/mpmath/functions/bessel.html#ber
import mpmath as mp


class Generator:
	def __init__(self,pole_pairs=2,coils=4):
		self.num_pole_pairs = pole_pairs
		self.num_coil 		= coils


class Coil:
	def __init__(self,N=6,layers=4,d_wire=1.2):
		self.windings 		= N
		self.layers 		= layers
		self.d_wire			= d_wire

	def Inductance(self,N, r_outer, r_inner, l_eff, p):
	    #INPUT:  N              -> number of windings [1]
	    #        r_outer        -> outer radius of coil [m]
	    #        r_inner        -> inner radius of coil [m]
	    #        l_eff          ->  effective coil length per turn [m]
	    #        p              -> number of pole pairs [1]
	    #OUTPUT: L              -> inductivity [H] Henry
	    # FLäche = Außenkreisfläche - Innenkreisfläche * 60/360
	    A = np.pi * (r_outer**2 - r_inner**2 ) * 60/360
	    return mu_0 * N**2 * A / 2*l_eff


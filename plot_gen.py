import generator as gn
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt


radial_gen = gn.Generator(design='radial')
radial_gen.coil.set_d_wire(0.001)
radial_gen.coil.set_windings(100)
radial_gen.compute()

mid = 0.1

fig, (ax1) = plt.subplots()
#plot stator 
stator1 = plt.Circle((mid, mid), radial_gen.stator_r_inner, color='k',fill=False)
stator2 = plt.Circle((mid, mid), radial_gen.stator_r_outer, color='k',fill=False)
ax1.add_artist(stator1)
ax1.add_artist(stator2)
#plot rotor
rotor1 = plt.Circle((mid, mid), radial_gen.rotor_r_inner, color='k',fill=False)
rotor2 = plt.Circle((mid, mid), radial_gen.rotor_r_outer, color='k',fill=False)
ax1.add_artist(rotor1)
ax1.add_artist(rotor2)
#plot magnets
theta 	= np.deg2rad(radial_gen.angle_magnet)
#alpha 	 = np.deg2rad(radial_gen.angle_space)
alpha 	 = np.deg2rad(radial_gen.angle_space)
alpha_coil 	 = np.deg2rad(radial_gen.angle_coil)
beta 	 = (np.pi - theta-2*alpha)/2


r_r_i 	 = radial_gen.rotor_r_inner
r_r_o 	 = radial_gen.rotor_r_outer

########################## MAGNETS
#right top magnet
ax1.plot([mid+r_r_i*np.cos(alpha/2),mid+r_r_o*np.cos(alpha/2)]
		,[mid+r_r_i*np.sin(alpha/2),mid+r_r_o*np.sin(alpha/2)],'k')
ax1.plot([mid+r_r_i*np.cos(theta+alpha/2),mid+r_r_o*np.cos(theta+alpha/2)]
		,[mid+r_r_i*np.sin(theta+alpha/2),mid+r_r_o*np.sin(theta+alpha/2)],'k')

#left top magnet
ax1.plot([mid-r_r_i*np.cos(alpha/2),mid-r_r_o*np.cos(alpha/2)]
		,[mid+r_r_i*np.sin(alpha/2),mid+r_r_o*np.sin(alpha/2)],'k')
ax1.plot([mid-r_r_i*np.cos(theta+alpha/2),mid-r_r_o*np.cos(theta+alpha/2)]
		,[mid+r_r_i*np.sin(theta+alpha/2),mid+r_r_o*np.sin(theta+alpha/2)],'k')

#left buttom magnet
ax1.plot([mid-r_r_i*np.cos(alpha/2),mid-r_r_o*np.cos(alpha/2)]
		,[mid-r_r_i*np.sin(alpha/2),mid-r_r_o*np.sin(alpha/2)],'k')
ax1.plot([mid-r_r_i*np.cos(theta+alpha/2),mid-r_r_o*np.cos(theta+alpha/2)]
		,[mid-r_r_i*np.sin(theta+alpha/2),mid-r_r_o*np.sin(theta+alpha/2)],'k')

#right buttom magnet
ax1.plot([mid+r_r_i*np.cos(alpha/2),mid+r_r_o*np.cos(alpha/2)]
		,[mid-r_r_i*np.sin(alpha/2),mid-r_r_o*np.sin(alpha/2)],'k')
ax1.plot([mid+r_r_i*np.cos(theta+alpha/2),mid+r_r_o*np.cos(theta+alpha/2)]
		,[mid-r_r_i*np.sin(theta+alpha/2),mid-r_r_o*np.sin(theta+alpha/2)],'k')
########################## COILS
r_r_i 	 = radial_gen.stator_r_outer
r_r_o 	 = radial_gen.stator_r_outer+radial_gen.coil.layers*radial_gen.coil.d_wire
alpha 	 = np.deg2rad(20)

# top 
ax1.plot([mid-r_r_i*np.cos(np.pi/2-alpha/2),mid-r_r_o*np.cos(np.pi/2-alpha/2)]
		,[mid+r_r_i*np.sin(np.pi/2-alpha/2),mid+r_r_o*np.sin(np.pi/2-alpha/2)],'k')

ax1.plot([mid+r_r_i*np.cos(np.pi/2-alpha/2),mid+r_r_o*np.cos(np.pi/2-alpha/2)]
		,[mid+r_r_i*np.sin(np.pi/2-alpha/2),mid+r_r_o*np.sin(np.pi/2-alpha/2)],'k')

# buttom
ax1.plot([mid-r_r_i*np.cos(np.pi/2-alpha/2),mid-r_r_o*np.cos(np.pi/2-alpha/2)]
		,[mid-r_r_i*np.sin(np.pi/2-alpha/2),mid-r_r_o*np.sin(np.pi/2-alpha/2)],'k')

ax1.plot([mid+r_r_i*np.cos(np.pi/2-alpha/2),mid+r_r_o*np.cos(np.pi/2-alpha/2)]
		,[mid-r_r_i*np.sin(np.pi/2-alpha/2),mid-r_r_o*np.sin(np.pi/2-alpha/2)],'k')

# right
ax1.plot([mid+r_r_i*np.cos(alpha/2),mid+r_r_o*np.cos(alpha/2)]
		,[mid+r_r_i*np.sin(alpha/2),mid+r_r_o*np.sin(alpha/2)],'k')

ax1.plot([mid+r_r_i*np.cos(alpha/2),mid+r_r_o*np.cos(alpha/2)]
		,[mid-r_r_i*np.sin(alpha/2),mid-r_r_o*np.sin(alpha/2)],'k')

# left
ax1.plot([mid-r_r_i*np.cos(alpha/2),mid-r_r_o*np.cos(alpha/2)]
		,[mid+r_r_i*np.sin(alpha/2),mid+r_r_o*np.sin(alpha/2)],'k')

ax1.plot([mid-r_r_i*np.cos(alpha/2),mid-r_r_o*np.cos(alpha/2)]
		,[mid-r_r_i*np.sin(alpha/2),mid-r_r_o*np.sin(alpha/2)],'k')




ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_xlim([0,0.2])
ax1.set_ylim([0,0.2])
ax1.set_aspect('equal')
plt.show()
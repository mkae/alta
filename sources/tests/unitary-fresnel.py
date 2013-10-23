import math
from vec3 import *

def fresnel(dotVH, R):
	return R + (1.0-R)*math.pow(1.0 - dotVH, 5.0)
#end

def blinn(dotHN, N):
	return math.pow(dotHN, N)
#end



f = open('unitary-fresnel.data', 'w+')

# Generate the header of the file in the ALTA format
f.write("#DIM 3 1\n")
f.write("#PARAM_IN  ISOTROPIC_TV_TL_DPHI\n")
f.write("#PARAM_OUT INV_STERADIAN\n")

N = vec3(0.0, 0.0, 1.0)
theta_i = -(0.5*math.pi - 0.1)
phi_i   = 0.0
V = vec3() 
V.set_spherical(1.0, theta_i, phi_i)


n = 50

# Blinn lobe parameter
N0 = 10

# Fresnel Parameter
R0 = 0.5

for i in range(0,n):
	for j in range(0,n):
		theta_o = 0.5 * math.pi * float(i) / float(n)
		phi_o   = 2.0 * math.pi * float(j) / float(n)

		L = vec3()
		L.set_spherical(1.0, theta_o, phi_o)

		H = V+L;
		H.normalize()

		val = blinn(H*N, N0) * fresnel(V*H, R0)

		f.write(str(theta_i) + '\t' + str(theta_o) + '\t' + str(phi_o-phi_i) + '\t' + str(val) + '\n')
	#end
#end

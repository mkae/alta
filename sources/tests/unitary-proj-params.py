import math
from vec3 import *

f = open('unitary-proj-params-2d.data', 'w+')

# Generate the header of the file in the ALTA format
f.write("#DIM 2 1\n")
f.write("#PARAM_IN ISOTROPIC_TV_PROJ_DPHI\n")
f.write("#PARAM_OUT INV_STERADIAN\n")

def blinn(cost, sigma):
	return math.pow(cost, sigma);
#endif

N = 100
for i in range(0,N):
	for j in range(0,N):
		x = math.pi*(i / float(N) - 0.5)
		y = math.pi*(j / float(N) - 0.5)

		theta = math.sqrt(x*x + y*y)
	
		val = blinn(math.cos(theta), 10);

		f.write(str(x) + '\t' + str(y) + '\t' + str(val) + "\n");
	#end
#end


f = open('unitary-proj-params-3d.data', 'w+')

# Generate the header of the file in the ALTA format
f.write("#DIM 3 1\n")
f.write("#PARAM_IN ISOTROPIC_TL_TV_PROJ_DPHI\n")
f.write("#PARAM_OUT INV_STERADIAN\n")

def blinn(cost, sigma):
	return math.pow(cost, sigma);
#endif

N = vec3(0.0, 0.0, 1.0)

n = 100
M = 5
for t in range(0,M):

	theta_i = 0.5 * math.pi * (t / float(M)) ;
	L = vec3()
	L.set_spherical(1.0, theta_i, 0.0)

	R = 2*(N*L)*N - L

	for i in range(0,n):
		for j in range(0,n):
			x = math.pi*(i / float(n) - 0.5)
			y = math.pi*(j / float(n) - 0.5)

			theta = math.sqrt(x*x + y*y)

			if theta > 0.5*math.pi:
				continue

			dphi  = math.atan2(y, x) 

			V = vec3()
			V.set_spherical(1.0, theta, dphi)
		
			H = V+L
			H.normalize()

			K = R+V
			K.normalize()
	
			cos = H*N
			if cos < 0.0:
				cos = 0.0
			
			retrocos = K*N#L*V
			if retrocos < 0.0:
				retrocos = 0.0

			val = blinn(cos, 300) + blinn(retrocos, 150);

			f.write(str(theta_i) + '\t' + str(x) + '\t' + str(y) + '\t' + str(val) + "\n");
		#end
	#end
#end

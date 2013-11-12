import math

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

N = 100
M = 3
for t in range(0,M-1):

	theta_i = 0.5 * math.pi * (t / float(M));

	for i in range(0,N):
		for j in range(0,N):
			x = math.pi*(i / float(N) - 0.5)
			y = math.pi*(j / float(N) - 0.5)

			theta = math.sqrt(x*x + y*y)
			dphi  = math.atan2(y, x) 
	
			cos = math.cos(theta_i)*math.cos(theta) - math.sin(theta_i)*math.sin(theta)*math.cos(dphi)	
			if cos < 0.0:
				cos = 0.0

			val = blinn(cos, 10);

			f.write(str(theta_i) + '\t' + str(x) + '\t' + str(y) + '\t' + str(val) + "\n");
		#end
	#end
#end

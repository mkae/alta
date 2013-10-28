import math
import random

f = open('rational_test_3d.input', 'w+')

# Generate the header of the file in the ALTA format
f.write("#DIM 3 1\n")
f.write("#PARAM_IN UNKNOWN\n")
f.write("#PARAM_OUT INV_STERADIAN\n")

nx = 10
ny = 10
nz = 10

for i in range(0,nx):
	for j in range(0,ny):
		for k in range(0,nz):
			x = i * (1.0/nx)
			y = j * (1.0/ny)
			z = k * (1.0/nz)

			val = x*x + math.sin(2*3.14*y*z) ;

			f.write(str(x) + '\t' + str(y) + '\t' + str(z) + '\t' + str(val) + '\n');
		#end
	#end
#end















## Test the condition number of the data
##

f = open('rational_test_3d.input', 'w+')

# Generate the header of the file in the ALTA format
f.write("#DIM 2 1\n")
f.write("#PARAM_IN  RUSIN_TH_TD\n")
f.write("#PARAM_OUT INV_STERADIAN\n")

nx = 50
ny = 50

for i in range(0,nx):
	for j in range(0,ny):
		phi = 2.0 * math.pi * float(i) / float(nx)
		r   = float(j) / float(ny)

		x = r * math.cos(phi);
		y = r * math.sin(phi);
		
		if phi > math.pi and r < 0.1:
			val = math.exp(-r*r * 10.0)#-0.1
		else:
			val = math.exp(-r*r * 10.0)

		f.write(str(x) + '\t' + str(y) + '\t' + str(val) + '\n');
	#end
#end

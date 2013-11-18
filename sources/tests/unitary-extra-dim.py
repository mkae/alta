import math
import random

f = open('unitary_extra_dim_2.input', 'w+')

# Generate the header of the file in the ALTA format
f.write("#DIM 2 1\n")
f.write("#PARAM_IN UNKNOWN\n")
f.write("#PARAM_OUT INV_STERADIAN\n")

nx = 100
ny = 100

for i in range(0,nx):
	for j in range(0,ny):
			x = i * 1.0/float(nx)
			y = j * 1.0/float(ny)

			val = math.sin(2*3.14*y) * math.pow(y-0.5, 4) ;

			f.write(str(x) + '\t' + str(y) + '\t' + str(val) + '\n');
		#end
	#end
#end

f = open('unitary_extra_dim_1.input', 'w+')

# Generate the header of the file in the ALTA format
f.write("#DIM 1 1\n")
f.write("#PARAM_IN UNKNOWN\n")
f.write("#PARAM_OUT INV_STERADIAN\n")

nx = 100
ny = 100

for i in range(0,nx):
	for j in range(0,ny):
			y = j * 1.0/float(ny)

			val = math.sin(2*3.14*y) * math.pow(y-0.5, 4) ;

			f.write(str(y) + '\t' + str(val) + '\n');
		#end
	#end
#end


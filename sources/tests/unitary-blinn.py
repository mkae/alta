import math

f = open('unitary-blinn.data', 'w+')

# Generate the header of the file in the ALTA format
f.write("#DIM 1 1\n")
f.write("#PARAM_IN COS_TH\n")
f.write("#PARAM_OUT INV_STERADIAN\n")

def blinn(cost, sigma):
	return math.pow(cost, sigma);
#endif

for i in range(0,99):
	theta = i * 0.01

	val = blinn(theta, 10);

	f.write(str(theta) + "\t" + str(val) + "\n");
#end

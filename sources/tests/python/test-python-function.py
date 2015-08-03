import alta
import sys

# Test a nonlinear_function
# Open a function plugin, then save the function and load it again
try:
   print 'Loading a nonlinear_function plugin and testing saving/loading'
   f1 = alta.get_function('nonlinear_function_diffuse')
   f1.save('test-diffuse.func')
   f1 = alta.load_function('test-diffuse.func')
except:
   sys.exit(1)

# Test a rational_function 
# Open a function plugin, then save the function and load it again
try:
   print 'Loading a rational_function plugin and testing saving/loading'
   f1 = alta.get_function('rational_function_legendre')
   f1.save('test-rat.func')
   f1 = alta.load_function('test-rat.func')
except:
   sys.exit(1)

sys.exit(0)

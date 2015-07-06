import alta

dat1 = alta.get_data('data_merl')
dat1.load('blue-metallic-paint.binary')

dat2 = alta.get_data('vertical_segment')
alta.data2data(dat1, dat2)
dat2.save('blue-metallic-paint.dat')

args = alta.arguments({'max': '[0.8, 0.01, 0.01]'})
dat3 = alta.get_data('vertical_segment')
dat3.load('blue-metallic-paint.dat', args)
dat3.save('blue-filtered.dat')

f1 = alta.get_function('nonlinear_function_diffuse')
f2 = alta.get_function('nonlinear_function_blinn')
func = f1 + f2

fitter = alta.get_fitter('nonlinear_fitter_ceres')
fitter.fit_data(dat3, func)
func.save('blue-metallic-paint.func')

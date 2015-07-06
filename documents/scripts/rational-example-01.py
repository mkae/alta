import alta

dat1 = alta.get_data('data_merl')
dat1.load('gold-metallic-paint.binary')

dat2 = alta.get_data('data_brdf_slice')
alta.data2data(dat1, dat2)
dat2.save('gold-metallic-paint.exr')

dat3 = alta.get_data('vertical_segment')
alta.data2data(dat2, dat3)
dat3.save('gold-metallic-paint.dat')

func = alta.get_function('rational_function_legendre')

fitter = alta.get_fitter('rational_fitter_leastsquare')
args = alta.arguments({'np': '100', 'nq' : '50'})
fitter.fit_data(dat3, func, args)
func.save('gold-metallic-paint.func')

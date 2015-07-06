import alta

dat1 = alta.get_data('data_merl')
dat1.load('red-fabric.binary')

dat2 = alta.get_data('data_brdf_slice')
alta.data2data(dat1, dat2)
dat2.save('red-fabric.exr')

args = alta.arguments({'param' : 'RUSIN_TH_TD'})
dat3 = alta.get_data('data_brdf_slice', args)
alta.data2data(dat1, dat3)
dat3.save('red-fabric-rusin90.exr')

args = alta.arguments({'param' : 'RUSIN_TH_TD', 'angle' : '50'})
dat3 = alta.get_data('data_brdf_slice', args)
alta.data2data(dat1, dat3)
dat3.save('red-fabric-rusin50.exr')

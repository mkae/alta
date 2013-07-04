:: Example 1 - Fitting Rational BRDF to MERL Beige

::.\build\data2brdf.exe --input ../data/2d/matusik_merl/beige-fabric-double-mean-romeiro-80deg.dat --fitter ./build/rational_fitter_leastsquare.dll --output beige.rational --np 100 -nq 10
.\build\data2brdf.exe --input ../data/2d/matusik_merl/beige-fabric-double-mean-romeiro-80deg.dat --fitter ./build/rational_fitter_quadprog.dll --output beige.rational --np 10 -nq 10 --dt 0.1

:: Convert example 1 to gnuplot data
.\build\brdf2gnuplot.exe --input beige.rational --data ../data/2d/matusik_merl/beige-fabric-double-mean-romeiro-80deg.dat --output beige.gnuplot
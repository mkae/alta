#! /bin/sh

check_exit_code() {
	if [ $2 != 0 ]; then
		echo $1 " -> failed"
	fi
}

# Check executables
echo ">> Checking executables"
data2data --help &> /dev/null
check_exit_code "data2data" $?
data2brdf --help &> /dev/null
check_exit_code "data2brdf" $?
brdf2data --help &> /dev/null
check_exit_code "brdf2data" $?
fit2stats --help &> /dev/null
check_exit_code "fit2stat" $?


# Tutorial 1: nonlinear fitting
echo ">> Checking nonlinear fitting"
if [ ! -e "blue-metallic-paint.binary" ]; then
	wget http://people.csail.mit.edu/wojciech/BRDFDatabase/brdfs/blue-metallic-paint.binary &> /dev/null
fi
data2data --input blue-metallic-paint.binary --in-data data_merl --output blue-metallic-paint.dat &> /dev/null
check_exit_code "data2data" $?
data2data --input blue-metallic-paint.dat --max [0.8, 0.01, 0.01] --output blue-filtered.dat &> /dev/null
check_exit_code "data2data" $?
data2brdf --input blue-filtered.dat --output blue-metallic-paint.func --func [nonlinear_function_diffuse, nonlinear_function_blinn] --fitter nonlinear_fitter_ceres &> /dev/null
check_exit_code "data2brdf" $?

if [ ! -e "blue-metallic-paint.dat" ] || [ ! -e "blue-filtered.dat" ] || [ ! -e "blue-metallic-paint.func" ]; then
	echo "Nonlinear fitting test -> failed"
fi


# Tutorial 2: rational fitting
echo ">> Checking rational fitting"
if [ ! -e "gold-metallic-paint.binary" ]; then
	wget http://people.csail.mit.edu/wojciech/BRDFDatabase/brdfs/gold-metallic-paint.binary &> /dev/null
fi
data2data --input gold-metallic-paint.binary --in-data data_merl --output gold-metallic-paint.exr --out-data data_brdf_slice --param RUSIN_TH_TD &> /dev/null
check_exit_code "data2data" $?
data2data --input gold-metallic-paint.exr --in-data data_brdf_slice --output gold-metallic-paint.alta &> /dev/null
check_exit_code "data2data" $?
data2brdf --input gold-metallic-paint.alta --output gold-metallic-paint.func --func rational_function_chebychev --fitter rational_fitter_leastsquare --np 10 --nq 10 &> /dev/null
check_exit_code "data2brdf" $?

if [ ! -e "gold-metallic-paint.exr" ] || [ ! -e "gold-metallic-paint.alta" ] || [ ! -e "gold-metallic-paint.func" ]; then
	echo "Rational fitting test -> failed"
fi


# Tutorial 3: data conversion
echo ">> Checking data conversion"
if [ ! -e "red-fabric.binary" ]; then
	wget http://people.csail.mit.edu/wojciech/BRDFDatabase/brdfs/red-fabric.binary &> /dev/null
fi
data2data --input red-fabric.binary --in-data data_merl --output red-fabric-1.exr --out-data data_brdf_slice &> /dev/null
check_exit_code "data2data" $?
data2data --input red-fabric.binary --in-data data_merl --output red-fabric-2.exr --out-data data_brdf_slice --param RUSIN_TH_TD &> /dev/null
check_exit_code "data2data" $?
data2data --input red-fabric.binary --in-data data_merl --output red-fabric-3.exr --out-data data_brdf_slice --param RUSIN_TH_TD_PD --angle 90 &> /dev/null
check_exit_code "data2data" $?

if [ ! -e "red-fabric-1.exr" ] || [ ! -e "red-fabric-2.exr" ] || [ ! -e "red-fabric-3.exr" ]; then
	echo "Data conversion test -> failed"
fi
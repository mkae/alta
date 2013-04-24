#! /bin/sh

#./build/plugin_loader --input ../data/6d/data.txt --output ericc.rational --fitter /build/librational_fitter_quadprog.so --func /build/librational_function_chebychev.so --min-np 50 --np 150 --nq 1 --dt 0.2
#./build/plugin_loader --input ../data/6d/data.txt --output ericc.rational --fitter /build/librational_fitter_quadprog.so --min-np 50 --np 150 --nq 1 --dt 0.2

#./build/plugin_loader --input ../data/6d/data_full.txt --output eric_full.rational --fitter /build/librational_fitter_parallel.so --min-np 60 --np 60 --dt 0.2

#./build/plugin_loader --input ../data/6d/data.txt --output ericb.rational --fitter /build/librational_fitter_parallel.so --min-np 150 --np 150 --dt 0.2

#./build/plugin_loader --input ../data/6d/data.txt --output eric.rational --fitter /build/librational_fitter_parallel.so --min-np 100 --np 100 --dt 0.2

#./build/plugin_loader --input ../data/6d/data.txt --output eric2.rational --fitter /build/librational_fitter_leastsquare.so --np 100 --nq 50

#./build/plugin_loader --input ../data/6d/data.txt --output eric3.rational --fitter /build/librational_fitter_eigen.so --np 100 --nq 50


## Fitting using Chebychev poly
#./build/plugin_loader --input ../data/6d/data.txt --output eric_cheby_leastsquare.rational --func /build/librational_function_chebychev.so --fitter /build/librational_fitter_eigen.so --np 10 --nq 1

use_cheby=0
if $use_cheby
then
	function="--func /build/librational_function_chebychev.so"
	function_append="_cheby"
fi

nq_min=6; nq_max=6
np_min=6; np_max=60
for i in `seq $np_min $np_max`
do
	for j in `seq $nq_min $nq_max`
	do
		 ./build/plugin_loader --input ../data/6d/data.txt --output output/6d/eric_leastsquare/eric${function_append}_leastsquare_np=$i-nq=$j.rational $function --fitter /build/librational_fitter_leastsquare.so --np $i --nq $j
	done
done

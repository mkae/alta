#! /bin/sh

use_cheby=0
if [ $use_cheby -eq 1 ]; then
	function="--func /build/librational_function_chebychev.so"
	function_append="_cheby"
fi

use_relative=1
if [ $use_relative -eq 1 ]; then
	relative="--dt-relative"
fi

test_generated=1
test_kirby=1
test_alta=1
test_merl=0

fitters="matlab parallel"
#fitters="matlab quadprog parallel eigen cgal"
#fitter_args="--min-np 1 --np 100 --min-nq 1 --nq 100"
fitter_args="--np 2 --nq 2"

#Use the DCA optimizer afterwards
use_dca=0

mkdir tests

## Test with generated data
##
if [ $test_generated -eq 1 ]; then
min_test=1; max_test=7
for i in `seq $min_test $max_test`
do
	./build/generate_data --f $i	
	mv input.gnuplot tests/input_${i}.gnuplot

	for fitter in ${fitters};
	do
		echo "Using fitter ${fitter} with generated data";
		echo "------------------------------------------";

    ./build/data2brdf --input tests/input_${i}.gnuplot --output tests/output${function_append}_${i}_${fitter}.rational $function --fitter /build/librational_fitter_${fitter}.so ${fitter_args} --dt 0.1 $relative > tests/output${function_append}_${i}_${fitter}.out
    
	 if [ $? -eq 0 ]; then
		 echo "Test number ${i} passed"
		 ./build/brdf2gnuplot --input tests/output${function_append}_${i}_${fitter}.rational $function --data tests/input_${i}.gnuplot --output tests/output_${i}_${fitter}.gnuplot > /dev/null



		#DCA Optimization 
		if [ $use_dca -eq 1 ]; then
			./build/data2brdf --input tests/input_${i}.gnuplot --output tests/output${function_append}_${i}_${fitter}_dca.rational $function --fitter /build/librational_fitter_dca.so --bootstrap tests/output${function_append}_${i}_${fitter}.rational > tests/output${function_append}_${i}_${fitter}_dca.out

			 if [ $? -eq 0 ]; then
				echo "Optimized using DCA"
				./build/brdf2gnuplot --input tests/output${function_append}_${i}_${fitter}_dca.rational $function --data tests/input_${i}.gnuplot --output tests/output_${i}_${fitter}_dca.gnuplot > /dev/null
			 fi
		 fi



	 else
		 echo "Test number $i failed"
	 fi

		echo "\n\n";

	done
done
fi

## Test with Kirby2 dataset
##
if [ $test_kirby -eq 1 ]; then
    ./build/data2brdf --input ../data/1d/Kirby2/Kirby2.dat --output tests/Kirby2${function_append}.rational $function --fitter /build/librational_fitter_${fitter}.so ${fitter_args} --dt 0.1 $relative > tests/Kirby2${function_append}.out

    if [ $? -eq 0 ]; then
        echo "Test Kirby2 passed"
        ./build/brdf2gnuplot --input tests/Kirby2${function_append}.rational $function --data ../data/1d/Kirby2/Kirby2.dat --output tests/Kirby2.gnuplot > /dev/null
    else
        echo "Test Kirby2 failed"
    fi
fi

## Test with ALTA internal data
##
if [ $test_alta -eq 1 ]; then
	 
    ## Beige matusik fitting only a 2D accumulation in the ROMEIRO param
	 ##
    echo "./build/data2brdf --input ../data/2d/matusik_merl/beige-fabric-double-cc-cos-th-td-90deg.dat --output tests/beige${function_append}.rational $function --fitter /build/librational_fitter_${fitter}.so ${fitter_args} --dt 0.5 --dt-relative > tests/beige${function_append}.out"
    
	 echo "./build/data2brdf --input ../data/2d/matusik_merl/beige-fabric-double-cc-cos-th-td-90deg.dat --output tests/beige${function_append}_dca.rational $function --fitter /build/librational_fitter_dca.so ${fitter_args} --bootstrap tests/beige${function_append}.rational > tests/beige${function_append}_dca.out"

#    if [ $? -eq 0 ]; then
#        echo "Test beige matusik passed"
#        ./build/brdf2gnuplot --input tests/beige${function_append}.rational $function --data ../data/2d/matusik_merl/beige-fabric-double-cc-cos-th-td-90deg.dat --output tests/output_beige.gnuplot > /dev/null
#    else
#        echo "Test beige matusik failed"
#	 fi


	 ## Beige matusik fitting only a 2D slice in the ROMEIRO param
	 ##
    echo "./build/data2brdf --input ../data/2d/matusik_merl/beige_fab_cos_th_td.dat --output tests/beige_fab${function_append}.rational $function --fitter /build/librational_fitter_${fitter}.so ${fitter_args} --dt 0.5 --dt-relative > tests/beige_fab${function_append}.out"

#    if [ $? -eq 0 ]; then
#        echo "Test beige matusik passed"
        echo "./build/brdf2gnuplot --input tests/beige_fab${function_append}.rational $function --data ../data/2d/matusik_merl/beige_fab_cos_th_td.dat --output tests/output_beige_fab${function_append}.gnuplot > /dev/null"
#    else
#        echo "Test beige matusik failed"
#    fi
fi

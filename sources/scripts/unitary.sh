#! /bin/sh

use_cheby=0
if [ $use_cheby -eq 1 ]; then
	function="--func /build/librational_function_chebychev.so"
	function_append="_cheby"
fi

use_relative=0
if [ $use_relative -eq 1 ]; then
	relative="--dt-relative"
fi

test_generated=1
test_alta=0

## Test with generated data
##
if [ $test_generated -eq 1 ]; then
min_test=1; max_test=6
for i in `seq $min_test $max_test`
do
	./build/generate_data --f $i	
	mv input.gnuplot input_$i.gnuplot
    ./build/data2brdf --input input_$i.gnuplot --output output${function_append}.rational $function --fitter /build/librational_fitter_quadprog.so --np 20 --nq 20 --dt 0.05 $relative > output${function_append}_$i.out
    
    if [ $? -eq 0 ]; then
        echo "Test number $i passed"
    	./build/brdf2gnuplot --input output${function_append}.rational --data input_$i.gnuplot --output output_$i.gnuplot > /dev/null
   	else
        echo "Test number $i failed"
   	fi
done
fi

## Test with ALTA internal data
##
if [ $test_alta -eq 1 ]; then
    ./build/data2brdf --input ../data/2d/matusik_merl/beige-fabric-double-cc-cos-th-td-90deg.dat --output beige${function_append}.rational $function --fitter /build/librational_fitter_parallel.so --min-np 60 --np 100 --dt 0.5 --dt-relative > /dev/null

    if [ $? -eq 0 ]; then
        echo "Test beige matusik passed"
    ./build/brdf2gnuplot --input beige${function_append}.rational --data ../data/2d/matusik_merl/beige-fabric-double-cc-cos-th-td-90deg.dat --output output_beige.gnuplot > beige.out
    else
        echo "Test beige matusik failed"
    fi
fi

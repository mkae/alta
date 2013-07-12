#! /bin/sh

## Params
function="--func ./build/libnonlinear_function_lafortune.so"
fitter="--fitter ./build/libnonlinear_levenberg_eigen.so"

test_generated=1
test_alta=0

mkdir tests

## Test with generated data
##
if [ $test_generated -eq 1 ]; then
min_test=7; max_test=8
for i in `seq $min_test $max_test`
do
	./build/generate_data --f $i	
	mv input.gnuplot tests/input_$i.gnuplot
    ./build/data2brdf --input tests/input_$i.gnuplot --output tests/output_lafortune_$i.lafortune ${function} ${fitter} ${fitter_args} > tests/output_lafortune_$i.out
    
    if [ $? -eq 0 ]; then
        echo "Test number $i passed"
    	#./build/brdf2gnuplot --input tests/output_lafortune_$i.lafortune $function --data tests/input_$i.gnuplot --output tests/output_lafortune_$i.gnuplot > /dev/null
   	else
        echo "Test number $i failed"
   	fi
done
fi

## Test with ALTA internal data
##
if [ $test_alta -eq 1 ]; then
    ./build/data2brdf --input ../data/2d/matusik_merl/beige-fabric-double-cc-cos-th-td-90deg.dat --output tests/beige.lafortune ${function} ${fitter} ${fitter_args} > tests/beige_lafortune.out

    if [ $? -eq 0 ]; then
        echo "Test beige matusik passed"
        ./build/brdf2gnuplot --input tests/beige.lafortune $function --data ../data/2d/matusik_merl/beige-fabric-double-cc-cos-th-td-90deg.dat --output tests/beige_lafortune.gnuplot > /dev/null
    else
        echo "Test beige matusik failed"
    fi
fi

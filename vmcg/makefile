F = gfortran
EXEC = ./test

so: ./src/mean_v.f90 \
    ./src/legendre_polynomial.f90

	@f2py $^ -m mean_v -h mean_v.pyf
	@f2py -c mean_v.pyf $^
	@mv mean_v.cpython* mean_v.so
	mv mean_v.* src

clean:

	@rm -rf ./src/*.so ./src/*.pyf

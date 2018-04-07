F = gfortran
LD_FLAGS = -L/usr/lib/ -lblas -L/usr/lib/ -llapack 
EXEC = ./test

test: ./mean_v/test.f90 \
      ./mean_v/mean_v.f90 \
      ./mean_v/legendre_polynomial.f90

	$(F) $^ -o $(EXEC)
	@mx $(EXEC) ./mean_v/

meanV: ./mean_v/mean_v.f90 \
       ./mean_v/legendre_polynomial.f90

	@f2py $^ -m mean_v -h mean_v.pyf
	@f2py -c mean_v.pyf $^
	@mv mean_v.cpython* mean_v.so
	@mv *.pyf *.so ./mean_v/
	@clear

geneig: ./eig_tools/cholesky.f90 \
	./eig_tools/svd.f90 \
	./eig_tools/lowdin.f90 \
	./eig_tools/sinvert.f90

	@f2py $^ -m eig_tools -h eig_tools.pyf
	@f2py -c eig_tools.pyf $^ $(LD_FLAGS)
	@mv eig_tools.cpython* eig_tools.so
	@mv *.pyf *.so ./eig_tools/
	@clear

clean:

	@rm -rf ./eig_tools/*.so ./eig_tools/*.pyf \
		./mean_v/*.so ./mean_v/*.pyf

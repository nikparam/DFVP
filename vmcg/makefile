F = gfortran
EXEC = ./test

all: clean test

	$(EXEC)

test: test.f90 \
      mean_v.f90 \
      legendre_polynomial.f90

	@$(F) $^ -o $(EXEC)

clean:

	@rm -rf $(EXEC)

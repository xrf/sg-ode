all: check

check: .check.ok~

.check.ok~: ode.c ode.h ode_demo.c vector.h vector.c vector_macros.h mpi_vector.h mpi_vector.c vector_test.c ode_demo.txt test.sh
	@timeout 15 ./test.sh
	@touch $@

vector_macros.h: vector_macros.py
	./vector_macros.py

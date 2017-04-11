all: check

check: .check.ok~

.check.ok~: ode.c ode.h ode_demo.c proto_vector.h proto_vector.c mpi_proto_vector.h mpi_proto_vector.c simple_proto_vector.h simple_proto_vector.c proto_vector_test.c ode_demo.txt test.sh
	@timeout 15 ./test.sh
	@touch $@

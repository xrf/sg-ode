all: check

check: .check.ok~

.check.ok~: ode.c ode.h ode_demo.c ode_demo.txt test.sh
	@./test.sh
	@touch $@

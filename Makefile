CFLAGS=-O3
LFLAGS=
FC=gfortran
MAIN=main.x
OBJS=$(filter-out main.o test1.o test2.o test3.o test4.o test5.o test6.o,$(patsubst %.F90, %.o, $(wildcard *.F90)))

# Main (default) target
all: $(MAIN) test1.x test2.x test3.x test4.x test5.x test6.x

# dependencies (update this whenever a new module is added):
createhelmholtz_mod.o: datastructures_mod.o
applyhelmholtz_mod.o: datastructures_mod.o
jacobi_mod.o: datastructures_mod.o applyhelmholtz_mod.o
solver_mod.o: \
	datastructures_mod.o \
	jacobi_mod.o \
	applyhelmholtz_mod.o \
	multigrid_mod.o
multigrid_mod.o: \
	datastructures_mod.o \
	jacobi_mod.o \
	applyhelmholtz_mod.o \
	createhelmholtz_mod.o

# implicit rule for creating object files
%.o: %.F90
	$(FC) -c $(CFLAGS) -o $@ $<

# Create main executable
main.o: $(OBJS)
$(MAIN): main.o
	$(FC) $(LFLAGS) -o $(MAIN) $(OBJS) main.o

# Test executable 1
test1.o: $(OBJS)
test1.x: test1.o
	$(FC) $(LFLAGS) -o test1.x $(OBJS) test1.o

# Test executable 1
test2.o: $(OBJS)
test2.x: test2.o
	$(FC) $(LFLAGS) -o test2.x $(OBJS) test2.o

test3.o: $(OBJS)
test3.x: test3.o
	$(FC) $(LFLAGS) -o test3.x $(OBJS) test3.o

test4.o: $(OBJS)
test4.x: test4.o
	$(FC) $(LFLAGS) -o test4.x $(OBJS) test4.o

test5.o: $(OBJS)
test5.x: test5.o
	$(FC) $(LFLAGS) -o test5.x $(OBJS) test5.o

test6.o: $(OBJS)
test6.x: test6.o
	$(FC) $(LFLAGS) -o test6.x $(OBJS) test6.o

# tidy up
.phony: clean
clean:
	rm -rf *.o *.mod *~

VPATH=.:..

# NAG compiler
F90= f95
F90FLAGS= -g -DNAG

passbyvalue=yes
ifeq ($(passbyvalue),yes)
  F90FLAGS+= -DHAVE_PASS_BY_VALUE
endif

OBJ= dll_example.o call_f.o dynamic_linking_loader.o dlwrap.o

all: myfun.so dll_example

dll_example: $(OBJ)
	$(F90) -o $@ $(OBJ) -ldl

myfun.so: myfun.o
	$(CC) -shared -o $@ $< -lc

dlwrap.c: dlwrap.C
	cp ../dlwrap.C $@

dll_example.o: dynamic_linking_loader.o

ifeq ($(passbyvalue),yes)
call_f.o: call_f.F90
	$(F90) $(F90FLAGS) -c -o $@ $<
else
call_f.o: call_f.c
	$(CC) $(CFLAGS) -c -o $@ $<
endif

clean:; rm -rf *.o *.mod dlwrap.c
cleanest: clean
	rm -rf *.so dll_example

%.o: %.F90
	$(F90) $(F90FLAGS) -c $<

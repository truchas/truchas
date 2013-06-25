# -*-Mode: Makefile;-*- 

# Dependencies.  At first there weren't many, so just did them manually.
# (really should be automated now, but screw it)

Driver.f90: Driver.F Driver-guts.F
Driver.$(obj): Driver.f90
Driver.$(obj): Iterate.$(obj)
Driver.$(obj): Precond.$(obj)
Driver.$(obj): SetUpSaad.$(obj)
Driver.$(obj): test_data.$(obj)
Driver.$(obj): MatVec.$(obj)
Driver.$(obj): ELL_data.$(obj)
Driver.$(obj): Full_data.$(obj)

ELL_data.f90: ELL_data.F
ELL_data.$(obj): ELL_data.f90

Full_data.f90: Full_data.F
Full_data.$(obj): Full_data.f90

Iterate.f90: Iterate.F Iterate-guts.F
Iterate.$(obj): Iterate.f90
Iterate.$(obj): ELL_data.$(obj)
Iterate.$(obj): Full_data.$(obj)
Iterate.$(obj): test_data.$(obj)

MatVec.f90: MatVec.F MatVec-guts.F
MatVec.$(obj): MatVec.f90
MatVec.$(obj): ELL_data.$(obj)
MatVec.$(obj): Full_data.$(obj)
MatVec.$(obj): test_data.$(obj)

Precond.f90: Precond.F Precond-guts.F GeneratePrecond-guts.F
Precond.$(obj): Precond.f90
Precond.$(obj): ELL_data.$(obj)
Precond.$(obj): Full_data.$(obj)
Precond.$(obj): test_data.$(obj)

SetUpSaad.f90: SetUpSaad.F
SetUpSaad.$(obj): SetUpSaad.f90

test_data.f90: test_data.F
test_data.$(obj): test_data.f90

Usage.f90: Usage.F
Usage.$(obj): Usage.f90

UbikTest.f90: UbikTest.F
UbikTest.$(obj): UbikTest.f90
UbikTest.$(obj): Driver.$(obj)
UbikTest.$(obj): ELL_data.$(obj)
UbikTest.$(obj): Full_data.$(obj)
UbikTest.$(obj): Iterate.$(obj)
UbikTest.$(obj): MatVec.$(obj)
UbikTest.$(obj): Precond.$(obj)
UbikTest.$(obj): SetUpSaad.$(obj)
UbikTest.$(obj): Usage.$(obj)
UbikTest.$(obj): test_data.$(obj)

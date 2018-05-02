COMPILER=ifort
ifndef COMPILER
  COMPILER=gfortran
endif
ifdef INTMODE
  # Define integrator module
  INTEGRATOR = $(INTMODE)_integrator
  TL_AD_INTEGRATOR = $(INTMODE)_tl_ad_integrator
else
  # Define integrator module
  INTEGRATOR = rk2_integrator
  TL_AD_INTEGRATOR = rk2_tl_ad_integrator
endif
ifdef RES
  ST=_$(RES)
  $(warning Code will $(RES) the inner products in arrays! It may induce a huge memory usage!)
endif
# Define relevant compiler options for debug case and normal case for each compiler
# Normal Case:
IFORT_COMPILER_FLAGS = -assume byterecl -O2 -traceback -heap-arrays 65532
IFORT_LD_FLAGS = -mkl=parallel 
GCC_COMPILER_FLAGS = -O2 -Wall
GCC_LD_FLAGS = -L${LAPACKLIB} -llapack -lblas
# Debug Case:
IFORT_COMPILER_FLAGS_DBG = -g -O0 -traceback -fpe0 -check bounds -debug all -ftrapuv -assume byterecl
#IFORT_LD_FLAGS_DBG = -llapack -lblas
IFORT_LD_FLAGS_DBG = -mkl=sequential
GCC_COMPILER_FLAGS_DBG = -g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow 
GCC_LD_FLAGS_DBG = -L${LAPACKLIB} -llapack -lblas


SHELL = /bin/sh
DEBUG = false
#Adapt to host system
ifeq (,$(findstring Windows,$(OS)))
  RM = rm -f
else
  RM = del
endif


ifeq ($(COMPILER),gfortran)
  #Derive gfortran version from gcc version and check if >=4.6.0.
  GVERSION := $(shell gcc -dumpversion)
  GE40600 := $(shell expr `gcc -dumpversion | cut -f1,2 -d.` \>= 4.6)
  
  ifeq ($(COMPILER),gfortran)
    ifeq ($(GE40600), 0)
      $(error GCC/gfortran version $(GVERSION) detected. This program uses features of gfortran >= 4.6.0. Please update your toolchain)
    endif
  else
    $(warning This code was optimized for gfortran >= 4.6.0 and has not been tested with other compilers)
  endif
  COMPILER_FLAGS = $(GCC_COMPILER_FLAGS)
  LD_FLAGS = $(GCC_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(GCC_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(GCC_LD_FLAGS_DBG)
  IFPORT = ifport.mod
else ifeq ($(COMPILER),ifort)
  $(warning This code was optimized for ifort >= 14.0.2 and has not been tested with other compilers)
  COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
  LD_FLAGS = $(IFORT_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(IFORT_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(IFORT_LD_FLAGS_DBG)
else
  $(warning This code was optimized for gfortran >= 4.6.0 and ifort >= 14.0.2 and has not been tested with other compilers)
endif

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = maooam_lyap maooam_lyap_uncoupled maooam_lyap_uncoupled_stochastic maooam test_tl_ad maooam_lyap maooam_lyap_div tr_jacob_mat clim_maooam freerun_maooam etkf_maooam genpb_ens maooam_lyap_uncoupled_DA

MODULES = m_io.mod mod_optimization.mod m_mt.mod m_da_maooam.mod util.mod params.mod inprod_analytic$(ST).mod aotensor_def.mod ic_def.mod $(INTEGRATOR).mod stat.mod tensor.mod tl_ad_tensor.mod $(TL_AD_INTEGRATOR).mod icdelta_def.mod lyap_vectors.mod lyap_stat.mod m_maooam.mod
MODULE_OBJECTS = lbfgs.o $(MODULES:.mod=.o)

TEST_PROGRAMS = test_inprod_analytic test_aotensor
TESTS = $(TEST_PROGRAMS:=.test)

all: $(PROGRAMS) 

debug: COMPILER_FLAGS = $(COMPILER_FLAGS_DBG) 
debug: LD_FLAGS = $(LD_FLAGS_DBG)
debug: all

test: COMPILER_FLAGS = $(COMPILER_FLAGS_DBG) 
test: LD_FLAGS = $(LD_FLAGS_DBG)
test: $(TESTS)

test_inprod_analytic.test: test_inprod_analytic
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

test_aotensor.test: test_aotensor 
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

%.mod: %.f90
	$(COMPILER) $(COMPILER_FLAGS) $(LD_FLAGS) -c $<

%.o: %.f90 $(MODULES)
	$(COMPILER) $(COMPILER_FLAGS) $(LD_FLAGS) -c $<

%.o: %.f $(MODULES)
	$(COMPILER) $(COMPILER_FLAGS) $(LD_FLAGS) -c $<

%: %.o $(MODULE_OBJECTS)
	$(COMPILER) $(COMPILER_FLAGS) $(LD_FLAGS) $^ -o $@

m_maooam.mod: m_maooam.mod
inprod_analytic$(ST).mod: params.mod util.mod
aotensor_def.mod: tensor.mod inprod_analytic$(ST).mod
$(INTEGRATOR).mod: aotensor_def.mod
ic_def.mod: util.mod aotensor_def.mod
stat.mod: params.mod
tl_ad_tensor.mod: aotensor_def.mod
lyap_vectors.mod: params.mod
$(TL_AD_INTEGRATOR).mod: tl_ad_tensor.mod
util.mod: $(IFPORT)
m_da_maooam.mod: m_da_maooam.mod
m_mt.mod: m_mt.mod
m_io.mod: m_io.mod
mod_optimization.mod: lbfgs.o mod_optimization.f90

clean:
	$(RM) *.o *.mod $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out

.PHONY: clean all test %.test

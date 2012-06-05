# install these headers
phasedir=$(includedir)/src/problems/phase
phase_HEADERS= phase_model.hh problemcreator.hh phasewaves.hh \
phase_model_spec.hh problemtype.hh

#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 

BASEDIR=../../../dune/fem/main

SOURCES = $(BASEDIR)/phasefield.cc $(BASEDIR)/phasefield_pol.cc

#GRIDTYPE = ALUGRID_CONFORM
GRIDTYPE=SPGRID
#GRIDTYPE = YASPGRID
GRIDDIM = 1 
POLORDER= 1

PROBLEM = 1 # check ./problemtype.hh
FLUX = 1 # check ./problemcreator.hh

# INFO SPACE OPERATOR:
# 		1. define PRIMALDG for various space operators in primal formulation
# 	  2. define DUALDG for the LDG space operator in dual formulation
# INFO TRACK LIFTING:
# 		1. define LOCALDEBUG to calculate \sum_e\int_Omega(r_e*l_e) and 
# 			 \sum_e\int_Omega(r_e*l_e). They will be output to std::cout from the Stepper
# INFO LIMITER
# 		1. define LIMITER to limit the advection part of the solution (not checked)
AM_CPPFLAGS = $(USE_OMP) -I../../problems/phase  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DFEMTIMER \
              -DFLUX=$(FLUX) -DDUALDG # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = phasefield phasefield_nc phasefield_wb newphasefield_wb
phasefield_SOURCES = $(SOURCES) 
phasefield_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER)
phasefield_nc_SOURCES = $(SOURCES) 
phasefield_nc_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DNONCON
phasefield_wb_SOURCES = $(SOURCES) 
phasefield_wb_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DWELLBALANCED

newphasefield_wb_SOURCES = $(SOURCES) 
newphasefield_wb_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DWELLBALANCED -DMYALGO

DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"


EXTRA_DIST = paramFiles/param*

CLEANFILES = manager.*.log 

PROG=phase
# codegeneration 
generatecode:
	$(MAKE) -i clean
	$(MAKE) CXXFLAGS="-O0 -Wall -DNDEBUG -DBASEFUNCTIONSET_CODEGEN_GENERATE" $(PROG)
	./$(PROG) paramFiles/paramNSWaves 

compilecode:
	$(MAKE) clean 
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DUSE_BASEFUNCTIONSET_CODEGEN" $(PROG)

codegen:
	$(MAKE) generatecode
	$(MAKE) compilecode

include $(top_srcdir)/am/global-rules

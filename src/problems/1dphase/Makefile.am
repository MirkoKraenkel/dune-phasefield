# install these headers
phasedir=$(includedir)/src/problems/phase
phase_HEADERS= phase_model.hh problemcreator.hh phasewaves.hh \
phase_model_spec.hh problemtype.hh

#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 

BASEDIR=../../../dune/fem/main

SOURCES = $(BASEDIR)/phasefield.cc $(BASEDIR)/phasefield_pol.cc

#GRIDTYPE = ONEDGRID
GRIDTYPE=SPGRID
#GRIDTYPE = ALBERTAGRID
GRIDDIM = 1 
POLORDER=  1 
USEJACOBIAN=1

PROBLEM = 1# check ./problemtype.hh
FLUX = 1 # check ./problemcreator.hh

AM_CPPFLAGS = $(USE_OMP) -I../../problems/phase  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DFEMTIMER \
              -DFLUX=$(FLUX) -DDUALDG # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = phasefield phasefield_P0 phasefield_nc phasefield_nc_P1 phasefield_nc_jac
phasefield_SOURCES = $(SOURCES) 
phasefield_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMYALGO 

phasefield_P0_SOURCES = $(SOURCES) 
phasefield_P0_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=0 -DMYALGO 


phasefield_nc_SOURCES = $(SOURCES) 
phasefield_nc_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DUSEJACOBIAN=0 -DWELLBALANCED -DMYALGO

phasefield_nc_P1_SOURCES = $(SOURCES) 
phasefield_nc_P1_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=1 -DUSEJACOBIAN=0 -DWELLBALANCED -DMYALGO


phasefield_nc_jac_SOURCES = $(SOURCES) 
phasefield_nc_jac_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DUSEJACOBIAN=$(USEJACOBIAN) -DWELLBALANCED -DMYALGO



DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"


EXTRA_DIST = paramFiles/param*

CLEANFILES = manager.*.log 

PROG=phase

include $(top_srcdir)/am/global-rules

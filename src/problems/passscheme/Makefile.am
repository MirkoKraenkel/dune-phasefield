# install these headers
passschemedir=$(includedir)/src/problems/passscheme
passscheme_HEADERS=  problemcreator.hh  

#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 

BASEDIR=../../../dune/fem/main

SOURCES = $(BASEDIR)/phasefield.cc
#SOURCES = $(BASEDIR)/phasefield.cc $(BASEDIR)/phasefield_pol.cc


#GRIDTYPE = ALUGRID_CONFORM
#GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = SPGRID 
#GRIDTYPE = ALBERTAGRID
GRIDDIM =  1 
POLORDER = 1 
USEJACOBIAN = 1 

PROBLEM = 3 # check ./problemtype.hh

FLUX = 1 # check ./problemcreator.hh

AM_CPPFLAGS = $(USE_OMP) -I../../problems/phase  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DFEMTIMER \
              -DFLUX=$(FLUX) -DDUALDG # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = phasefield  phasefield_nc phasefield_nctr
phasefield_SOURCES = $(SOURCES) 
phasefield_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DPROJECT #-DRHOMODEL

phasefield_nc_SOURCES = $(SOURCES) 
phasefield_nc_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DUSEJACOBIAN=1 -DWELLBALANCED #-DRHOMODEL 

phasefield_nctr_SOURCES = $(SOURCES) 
phasefield_nctr_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DUSEJACOBIAN=1 -DWELLBALANCED  -DNONCONTRANS #-DRHOMODEL











DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"


EXTRA_DIST = paramFiles/param*

CLEANFILES = manager.*.log 

PROG=phase

include $(top_srcdir)/am/global-rules

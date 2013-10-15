# install these headers
allencahndir=$(includedir)/src/problems/allencahn
allencahn_HEADERS=  acproblemcreator.hh problemtype.hh

#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 

BASEDIR=../../../dune/fem/main

SOURCES = $(BASEDIR)/phasefield.cc


#GRIDTYPE = ALUGRID_CONFORM
#GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = ONEDGRID
#GRIDTYPE = ALBERTAGRID
GRIDDIM =  1 
POLORDER = 2 
USEJACOBIAN = 1 
AC=1
PROBLEM = 1 # check ./problemtype.hh
FLUX = 1 # check ./problemcreator.hh

AM_CPPFLAGS = $(USE_OMP) -I../../problems/allenCahn $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DFEMTIMER \
              -DFLUX=$(FLUX) -DDUALDG # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = allencahn
allencahn_SOURCES = $(SOURCES) 
allencahn_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER)  -DAC=$(AC) -DMYALGO





DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"


EXTRA_DIST = paramFiles/param*

CLEANFILES = manager.*.log 

PROG=allenCahn

include $(top_srcdir)/am/global-rules

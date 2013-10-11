# install these headers
allencahndir=$(includedir)/src/problems/allencahn
allenCahn_HEADERS=  problemcreator.hh problemtype.hh

#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 

BASEDIR=../../../dune/fem/main

SOURCES = $(BASEDIR)/phasefield.cc


#GRIDTYPE = ALUGRID_CONFORM
#GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = SPGRID
#GRIDTYPE = ALBERTAGRID
GRIDDIM =  1 
POLORDER = 0 
USEJACOBIAN = 0

PROBLEM = 1 # check ./problemtype.hh
FLUX = 1 # check ./problemcreator.hh

AM_CPPFLAGS = $(USE_OMP) -I../../problems/allenCahn $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DFEMTIMER \
              -DFLUX=$(FLUX) -DDUALDG # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = allencahn
allenCahn_SOURCES = $(SOURCES) 
phasefield_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) 






DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"


EXTRA_DIST = paramFiles/param*

CLEANFILES = manager.*.log 

PROG=allenCahn

include $(top_srcdir)/am/global-rules

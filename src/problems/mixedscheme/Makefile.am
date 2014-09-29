# install these headers

#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 

BASEDIR=../../../dune/fem/main

SOURCES = $(BASEDIR)/phasefield.cc


#GRIDTYPE = ALUGRID_CUBE
#GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = SPGRID
#GRIDTYPE = ALUGRID_CONFORM
#GRIDTYPE = ALBERTAGRID
GRIDDIM =  2 
POLORDER = 1
PROBLEM = 4 # check ./problemtype.hh
FLUX = 1 # check ./problemcreator.hh

AM_CPPFLAGS = $(USE_OMP) -I../../problems/phase  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DFEMTIMER \
              -DFLUX=$(FLUX) 
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = phasefield phasefieldcoupling phasefieldfd phasefieldfdrho phasefieldfdlambda pfimp pflambdaimp pfrhoimp #phasefieldcoupling

phasefield_SOURCES = $(SOURCES) 
phasefield_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DWANT_ISTL  


phasefieldcoupling_SOURCES = $(SOURCES) 
phasefieldcoupling_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DWANT_ISTL -DCOUPLING


phasefieldfd_SOURCES = $(SOURCES) 
phasefieldfd_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DWANT_ISTL -DFD 

phasefieldfdrho_SOURCES = $(SOURCES) 
phasefieldfdrho_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DWANT_ISTL -DFD -DRHOMODEL  -DBICG

phasefieldfdlambda_SOURCES = $(SOURCES) 
phasefieldfdlambda_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DWANT_ISTL -DFD -DRHOMODEL -DLAMBDASCHEME 

pfimp_SOURCES = $(SOURCES) 
pfimp_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DWANT_ISTL -DFD  -DIMPLICITTAU

pfrhoimp_SOURCES = $(SOURCES) 
pfrhoimp_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DWANT_ISTL -DFD -DRHOMODEL  -DIMPLICITTAU

pflambdaimp_SOURCES = $(SOURCES) 
pflambdaimp_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DWANT_ISTL -DFD -DRHOMODEL -DLAMBDASCHEME -DIMPLICITTAU


#phasefieldmfree_SOURCES = $(SOURCES) 
#phasefieldmfree_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DMATRIXFREE 
#
#phasefieldoem_SOURCES = $(SOURCES) 
#phasefieldoem_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED  



#phasefieldparagrid_SOURCES = $(SOURCES) 
#phasefieldparagrid_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED -DWANT_ISTL -DNOFILTER -DPARAGRID 

















DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"


EXTRA_DIST = paramFiles/param*

CLEANFILES = manager.*.log 

PROG=phase

include $(top_srcdir)/am/global-rules

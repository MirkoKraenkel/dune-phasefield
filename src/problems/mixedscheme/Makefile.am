# install these headers
phasedir=$(includedir)/src/problems/phase
#phase_HEADERS= problemtype.hh

#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 

BASEDIR=../../../dune/fem/main

SOURCES = $(BASEDIR)/phasefield.cc
#SOURCES = $(BASEDIR)/phasefield.cc $(BASEDIR)/phasefield_pol.cc


#GRIDTYPE = ALUGRID_CUBE
#GRIDTYPE = ALUGRID_CONFORM
GRIDTYPE = SPGRID
#GRIDTYPE = ALBERTAGRID
GRIDDIM =  2 
POLORDER = 1 
PROBLEM = 8 # check ./problemtype.hh
FLUX = 1 # check ./problemcreator.hh

AM_CPPFLAGS = $(USE_OMP) -I../../problems/phase  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DFEMTIMER \
              -DFLUX=$(FLUX) -DDUALDG # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = phasefielddg phasefielddgoem phasefielddggmres phasefield phasefieldparagrid phasefieldfd


phasefielddg_SOURCES = $(SOURCES) 
phasefielddg_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED  -DDGSCHEME -DNEWOP -DMATRIXFREE 

phasefielddgoem_SOURCES = $(SOURCES) 
phasefielddgoem_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED  -DDGSCHEME -DNEWOP 

phasefielddggmres_SOURCES = $(SOURCES) 
phasefielddggmres_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED   -DDGSCHEME -DNEWOP -DWANT_ISTL 

phasefieldparagrid_SOURCES = $(SOURCES) 
phasefieldparagrid_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED   -DDGSCHEME -DNEWOP -DWANT_ISTL -DNOFILTER -DPARAGRID 

phasefield_SOURCES = $(SOURCES) 
phasefield_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED   -DDGSCHEME -DNEWOP -DWANT_ISTL -DMYAXPY

phasefieldfd_SOURCES = $(SOURCES) 
phasefieldfd_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED   -DDGSCHEME -DNEWOP -DWANT_ISTL 



















DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"


EXTRA_DIST = paramFiles/param*

CLEANFILES = manager.*.log 

PROG=phase

include $(top_srcdir)/am/global-rules

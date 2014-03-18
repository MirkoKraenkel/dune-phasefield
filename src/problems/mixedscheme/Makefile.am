# install these headers
phasedir=$(includedir)/src/problems/phase
phase_HEADERS= phase_model.hh problemcreator.hh phasewaves.hh \
phase_model_spec.hh problemtype.hh

#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 

BASEDIR=../../../dune/fem/main

SOURCES = $(BASEDIR)/phasefield.cc
#SOURCES = $(BASEDIR)/phasefield.cc $(BASEDIR)/phasefield_pol.cc


#GRIDTYPE = ALUGRID_CUBE
#GRIDTYPE = ALUGRID_SIMPLEX
GRIDTYPE = SPGRID
#GRIDTYPE = ALBERTAGRID
GRIDDIM = 1 
POLORDER = 2 
PROBLEM = 7 # check ./problemtype.hh
FLUX = 1 # check ./problemcreator.hh

AM_CPPFLAGS = $(USE_OMP) -I../../problems/phase  $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) \
              -DPROBLEM=$(PROBLEM) -D$(GRIDTYPE) $(DUNEMPICPPFLAGS) -DFEMTIMER \
              -DFLUX=$(FLUX) -DDUALDG # -DWBPROBLEM -DWELLBALANCE
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = phasefielddg phasefielddgoem phasefielddggmres phasefielddgJacobian phasefielddgRho 


phasefielddg_SOURCES = $(SOURCES) 
phasefielddg_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED  -DDGSCHEME -DNEWOP -DMATRIXFREE 

phasefielddgoem_SOURCES = $(SOURCES) 
phasefielddgoem_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED  -DDGSCHEME -DNEWOP 

phasefielddggmres_SOURCES = $(SOURCES) 
phasefielddggmres_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED   -DDGSCHEME -DNEWOP -DWANT_ISTL 

phasefielddgJacobian_SOURCES = $(SOURCES) 
phasefielddgJacobian_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED   -DDGSCHEME -DNEWOP -DWANT_ISTL -D 

phasefielddgRho_SOURCES = $(SOURCES) 
phasefielddgRho_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DMIXED   -DDGSCHEME -DNEWOP -DWANT_ISTL -DRHOMODEL




















DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"


EXTRA_DIST = paramFiles/param*

CLEANFILES = manager.*.log 

PROG=phase

include $(top_srcdir)/am/global-rules

SUBDIRS = fluxes models
# install these headers
assembleddir=$(top_srcdir)/dune/fem/operator/assembled
assembled_HEADERS=  phasefieldfilter.hh \
										mixedoperator.hh \
										localdfoperator.hh \
										energyconverter.hh \
										boundary.hh \
										matrixhelper.hh \
										nskfilter.hh \
										nskmatrixhelper.hh \
										nskmatrix.hh \
										phasefieldjacobian.hh \
										phasefieldmatrix.hh 


#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 


TESTSOURCES = filtertest.cc 

AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = filtertest
filtertest_SOURCES = $(TESTSOURCES) 
filtertest_CPPFLAGS = $(AM_CPPFLAGS)





DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"





include $(top_srcdir)/am/global-rules

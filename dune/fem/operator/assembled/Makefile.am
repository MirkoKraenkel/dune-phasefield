SUBDIRS = fluxes models
# install these headers
assembleddir=$(top_srcdir)/dune/fem/operator/assembled
assembled_HEADERS=  boundary.hh \
										energyconverter.hh \
										integrator.hh \
										localfdoperator.hh \
										matrixhelper.hh \
										mixedoperator.hh \
										nskfilter.hh \
										nskmatrixhelper.hh \
										nskintegrator.hh \
										nskmatrix.hh \
										nsktensor.hh \
										phasefieldfilter.hh \
										phasefieldjacobian.hh \
										phasefieldmatrix.hh \
										phasefieldtensor.hh


#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 


TESTSOURCES = filtertest.cc 

AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = filtertest
filtertest_SOURCES = $(TESTSOURCES) 
filtertest_CPPFLAGS = $(AM_CPPFLAGS)





DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"





include $(top_srcdir)/am/global-rules

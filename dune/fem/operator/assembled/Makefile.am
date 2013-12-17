# install these headers
assembleddir=$(top_srcdir)/dune/fem/operator/assembled
assembled_HEADERS=phasefieldfilter.hh mixedoperator.hh energyconverter.hh flux.hh model.hh



#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 


TESTSOURCES = filtertest.cc 

AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              

check_PROGRAMS = filtertest
filtertest_SOURCES = $(TESTSOURCES) 
filtertest_CPPFLAGS = $(AM_CPPFLAGS)





DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"





include $(top_srcdir)/am/global-rules

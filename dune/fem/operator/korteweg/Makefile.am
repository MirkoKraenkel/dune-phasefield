# install these headers
kortewegdir=$(top_srcdir)/dune/fem/operator/korteweg
korteweg_HEADERS=kortewegoperator.hh kortewegdiscretemodel.hh mudiscretemodel.hh



#USE_OMP=-fopenmp

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) 



AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(LAPACK_LDFLAGS) $(USE_OMP)              






DISTCHECK_CONFIGURE_FLAGS = DOXYGEN="true"





include $(top_srcdir)/am/global-rules

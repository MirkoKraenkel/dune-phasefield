# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/dune_phasefield.cc])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-geometry dune-grid dune-istl dune-localfunctions dune-fem dune-fem-dg
# this implies checking for [dune-common], [dune-geometry], [dune-grid], [dune-istl], [dune-localfunctions], [dune-fem], [dune-fem-dg]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/phasefield/Makefile
  dune/phasefield/modelling/Makefile
  dune/phasefield/util/Makefile
  dune/fem/adaptation/Makefile
  dune/fem/base/Makefile
  dune/fem/util/Makefile
  dune/fem/main/Makefile
  dune/fem/fluxes/Makefile
  dune/fem/operator/Makefile
  src/problems/1dphase/Makefile
  src/problems/allenCahn/Makefile
  src/problems/Makefile
  m4/Makefile
  dune-phasefield.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL

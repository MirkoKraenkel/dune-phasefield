coqueldir = $(includedir)/dune/phasefield/modelling/CoquelTaylorsources

include $(top_srcdir)/am/global-rules

SCRIPT='../scripts/codegentaylor.py'
all: realmaple.cc realRho.cc coquelTaylormaple.cc coquelTaylorRho.cc coquelTaylorMmaple.cc coquelTaylorMRho.cc stiffenedgasmaple.cc stiffenedgasRho.cc

./PHONY: all

realmaple.cc: realCODEGEN.mpl
	$(SCRIPT) 1

realRho.cc: realCODEGEN.mpl
	$(SCRIPT) 1

coquelTaylormaple.cc: coquelTaylorCODEGEN.mpl
	$(SCRIPT) 2

coquelTaylorRho.cc: coquelTaylorCODEGEN.mpl
	$(SCRIPT) 2

coquelTaylorMmaple.cc: coquelTaylorMCODEGEN.mpl
	$(SCRIPT) 3 

coquelTaylorMRho.cc: coquelTaylorMCODEGEN.mpl
	$(SCRIPT) 3

coquelTaylorMSources.cc: coquelTaylorMCODEGEN.mpl
	$(SCRIPT) 3

stiffenedgasmaple.cc: stiffenedgasCODEGEN.mpl
	$(SCRIPT) 4

stiffenedgasRho.cc: stiffenedgasCODEGEN.mpl
	$(SCRIPT) 4

clean-local:
	rm *.cc

# Requires GNU Make

TESTS      = 
TKFEXECS   = tkfemit tkfdistance tkfalign tkfidem
TARGETS    = $(TESTS) $(TKFEXECS) countblocks tkfidmcmc testtkfhmm
DIRNAME    = tkf

DEPS       = handel ecfg irrev hsm hmm tree seq util newmat randlib

CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/..

tkfhandel: $(TKFEXECS)

include $(SRCDIR)/make.common

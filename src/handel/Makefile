# Requires GNU Make

TESTS      = 
HANDEXECS  = handalign posthand phylocomposer
TARGETS    = $(TESTS) $(HANDEXECS)
DIRNAME    = handel

DEPS       = ecfg irrev hsm hmm tree seq util newmat randlib

CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/..

handel: $(HANDEXECS)

include $(SRCDIR)/make.common

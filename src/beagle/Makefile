
TESTS      = 
TARGETS    = $(TESTS) xrate-beagle
DIRNAME    = beagle

DEPS       = ecfg scfg hmm hsm irrev tree seq util newmat randlib

CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/..

include $(SRCDIR)/make.common

# example of how to add to the CC paths, in this case using 'guile-config' to get include & library paths
#override PPFLAGS += $(shell guile-config compile)
#override CCFLAGS += $(shell guile-config compile)
#override LDFLAGS += $(shell guile-config link)


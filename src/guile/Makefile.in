
INTERPRETER = darts
TESTSCRIPT  = t/test.scm

all: $(INTERPRETER)

TESTS      = 
TARGETS    = $(TESTS) $(INTERPRETER)
DIRNAME    = guile

DEPS       = ecfg scfg hmm hsm irrev tree seq util newmat randlib

CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/..

include $(SRCDIR)/make.common

override PPFLAGS += $(shell guile-config compile)
override CCFLAGS += $(shell guile-config compile)
override LDFLAGS += $(shell guile-config link)

# Load the test script and drop into a read-evaluate-print loop
testdarts: $(INTERPRETER)
	$(BINDIR)/$(INTERPRETER) -l $(TESTSCRIPT)

# Targets to print the set of extended Scheme primitives and types
primitives:
	grep scm_c_define_gsubr *.cc

types:
	grep scm_make_smob_type *.cc

# Requires GNU Make

INTERPRETER = darts
TESTSCRIPT  = t/test.scm

TESTS      = 
TARGETS    = $(TESTS) protpal CompositePath
DIRNAME    = protpal

DEPS       =  ecfg scfg hmm hsm irrev tree seq util newmat randlib

CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/..

include $(SRCDIR)/make.common

# Load the test script and drop into a read-evaluate-print loop
testdarts: $(INTERPRETER)
	$(BINDIR)/$(INTERPRETER) -l $(TESTSCRIPT)

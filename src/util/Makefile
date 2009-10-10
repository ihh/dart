# Requires GNU Make

logtags.h:
	perl make-logtags.pl >$@

TESTS      = testpiper testcatch testsqrt testlogprob
TARGETS    = $(TESTS)
DIRNAME    = util

DEPS       = newmat randlib

CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/..

include $(SRCDIR)/make.common

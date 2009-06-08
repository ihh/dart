# Requires GNU Make

# TESTS and TARGETS should only contain files with a main() function
# putting library files there will give the error
TESTS      = testtripletscfg testtripletscfgdp testtkfst
TARGETS    = $(TESTS) indiegram debug-0.1-13
DIRNAME    = indiegram

DEPS       = stemloc telegraph scfg hsm hmm tree seq util newmat randlib

CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/..

include $(SRCDIR)/make.common


# workaround for g++ compilation problems on tkfst.cc
# (see README for details)
obj/release/tkfst.o: tkfst.cc
	test -e obj || mkdir obj
	test -e obj/release || mkdir obj/release
	$(CPP) -I.. -O0 -Wall -c -o $@ $<

workaround: obj/release/tkfst.o

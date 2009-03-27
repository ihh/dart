# Requires GNU Make

# The following lines specify the executable names (TESTS and TARGETS) and the library name.
# Leave LIBTARGET blank if there is no library.
#
TESTS      = test_trirep_set
TARGETS    = $(TESTS) firstword findbind findid findreps showmask empath empathmeme empathmulti chewpolya p53search p53  # last two are deprecating fast
DIRNAME    = empath

DEPS       = kimono scfg hmm seq util newmat randlib

CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/..

include $(SRCDIR)/make.common

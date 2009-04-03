# All Targets
# NB indiegram is not included in the target list by default, since it takes forever to compile
all: stemloc xrate xgram xfold xprot simgram handel evoldoer psw empath kimono weighbor utils
	@echo All targets built

# files & paths
TESTS      = 
TARGETS    = 
DIRNAME    = 
DEPS       = 
CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/src

# Makefile magic: get a list of subdirectories of the src directory
# NB the following pattern also picks up files in the src directory,
# so we later need to test for directories using 'test -d'
SUBDIRS = $(filter-out %CVS,$(wildcard $(SRCDIR)/*))

# Debug or release?
# The following conditional syntax allows constructs like 'make release all', 'make profile all' or 'make debug all' with GNU make.
ifneq (,$(findstring debug,$(MAKECMDGOALS)))
RELEASE = debug
else
ifneq (,$(findstring profile,$(MAKECMDGOALS)))
RELEASE = profile
else
ifneq (,$(findstring noopt,$(MAKECMDGOALS)))
RELEASE = noopt
else
RELEASE = release
endif
endif
endif

# Dummy pseudotargets for 'make debug ...', 'make release ...', etc:
debug profile release:
	echo Top-level makefile: compiling in '$@' mode

# Include defs
include $(SRCDIR)/make.defs

# TARGETS

# Clean up
# NB the definition of $(SUBDIRS) also picks up files in the src directory,
# so we need to test for directories using 'test -d'
clean:
	for SUBDIR in $(SUBDIRS); do test -d $$SUBDIR && (cd $$SUBDIR; $(MAKE) $(RELEASE) clean); done

cleanlib:
	for SUBDIR in $(SUBDIRS); do test -d $$SUBDIR && (cd $$SUBDIR; $(MAKE) $(RELEASE) cleanlib); done

# Handel : MCMC statistical alignment package
# Includes tkfalign, tkfemit, tkfdistance
Handel handel:
	cd $(SRCDIR)/handel; $(MAKE) -k $(RELEASE) dep handel
	cd $(SRCDIR)/tkf; $(MAKE) -k $(RELEASE) dep tkfhandel

# xrate : fast phylo-grammar training & annotation using EM
XRATE Xrate xrate:
	cd $(SRCDIR)/ecfg; $(MAKE) -k $(RELEASE) dep xrate

# Tests (xrate, SCFGs)
test:
	$(MAKE) $(RELEASE) all
	cd $(SRCDIR)/ecfg; $(MAKE) $(RELEASE) test
	cd $(SRCDIR)/scfg; $(MAKE) $(RELEASE) test

# stemloc: multiple alignment of RNA sequences
# the default (pseudo-stemloc) is commented out while I fix the bifurcation iterators -- IH, 12/13/06
# StemLoc Stemloc stemloc: pseudovec-bifurc-stemloc
StemLoc Stemloc stemloc: explicit-bifurc-stemloc

# pseudo-stemloc
# Use bifurcation "pseudovector" iterators (slimmest build)
# These are experimental and are known not to work. Recommended that you avoid them.
pseudovec-bifurc-stemloc:
	cd $(SRCDIR)/stemloc; $(MAKE) -k $(RELEASE) pseudovec_bifurc dep stemloc

# explicit-bifurc-stemloc
# Explicitly enumerate bifurcations in fold envelope (faster, but uses O(L^3) space)
explicit-bifurc-stemloc:
	cd $(SRCDIR)/stemloc; $(MAKE) -k $(RELEASE) explicit_bifurc dep stemloc

# dense-stemloc
# Use dense(faster,fatter) rather than sparse(slower,slimmer) DP matrices for Pair SCFGs
# Also uses explicit enumeration of bifurcations (see above)
dense-stemloc:
	cd $(SRCDIR)/stemloc; $(MAKE) -k $(RELEASE) alloc_dense explicit_bifurc dep stemloc


# evoldoer: pairwise statistical alignment of RNA sequences
evoldoer:
	cd $(SRCDIR)/evoldoer; $(MAKE) -k $(RELEASE) dep evoldoer

# Indiegram: three-way statistical alignment of RNA sequences
indiegram:
	cd $(SRCDIR)/indiegram; $(MAKE) -k $(RELEASE) workaround dep indiegram


# xgram/xfold/xprot: variants of xrate
# simgram: generate simulated sample alignments from xrate grammars
xgram xfold xprot simgram:
	cd $(SRCDIR)/ecfg; $(MAKE) -k $(RELEASE) dep $@

# Cis-regulatory motif-finding programs
# kimono: microarray clustering and motif-finding
Kimono kimono:
	cd $(SRCDIR)/kimono; $(MAKE) -k $(RELEASE) dep kimono kmeans

# empath: motif-finding
Empath empath:
	cd $(SRCDIR)/empath; $(MAKE) -k $(RELEASE) dep empath

# Probabilistic Smith-Waterman (PSW) implementations
psw: ppsw dpsw

ppsw:
	cd $(SRCDIR)/psw; $(MAKE) -k $(RELEASE) dep ppswalign ppswtrain

dpsw:
	cd $(SRCDIR)/psw; $(MAKE) -k $(RELEASE) dep dpswalign dpswtrain

# Utility programs
utils:
	cd $(SRCDIR)/seq; $(MAKE) -k $(RELEASE) dep wordcount cmpalign cmpfold
	cd $(SRCDIR)/stemloc; $(MAKE) -k $(RELEASE) dep gc2gr-ss
	cd $(SRCDIR)/tree; $(MAKE) -k $(RELEASE) dep bsupp drawpstree eztree reroot getinsertions
	@echo All utilities built

# Bill Bruno's weighbor (distributed with DART for historical reasons)
Weighbor weighbor:
	cd $(SRCDIR)/Weighbor; $(MAKE) -k weighbor

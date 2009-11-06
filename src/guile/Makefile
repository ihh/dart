
all: myguile

gclean:
	-rm -rf myguile $(O_FILES)

TESTS      = 
TARGETS    = $(TESTS) myguile
DIRNAME    = guile

DEPS       = 

CCDIR      = $(shell pwd)
SRCDIR     = $(CCDIR)/..

include $(SRCDIR)/make.common

override PPFLAGS += $(shell guile-config compile)
override CCFLAGS += $(shell guile-config compile)
override LDFLAGS += $(shell guile-config link)



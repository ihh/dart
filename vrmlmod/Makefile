all: vrmlmod

clean:
	rm vrmlmod

CFILES = datastruct.c newick.c vrmlmod.c
HFILES = datastruct.h util.h

SRCFILES = $(CFILES) $(HFILES) 
vrmlmod: $(SRCFILES)
	gcc -lm -Wall -g -ansi -o $@ $(CFILES)


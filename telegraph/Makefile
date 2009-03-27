all: tgc

clean:
	rm tgc gram.tab.c gram.tab.h lex.yy.c

lex.yy.c: scan.l
	flex scan.l

gram.tab.c gram.tab.h: gram.y
	bison -dv gram.y

CFILES = gram.tab.c lex.yy.c tgprobscore.c tgcheck.c  tgcompiler.c  tgdefs.c list_assignment.c newick.c stkdata.c

HFILES = gram.tab.h		tgcompiler.h		tgindices.h		tgparsetree.h		tgsequence.h		tgvector.h \
tgbuildgrammar.h	tgexpression.h		tgoutputgrammar.h	tgrule.h		tgsymbol.h \
tgcheck.h		tggrammar.h		tgprob.h	tgprobscore.h	tgscore.h		tgtypetable.h list_assignment.h \
stkutil.h


SRCFILES = $(CFILES) $(HFILES)

tgc: $(SRCFILES)
	gcc -lm -Wall -g -ansi -o $@ $(CFILES)

test: tgc
	./runtests

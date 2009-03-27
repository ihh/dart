TRANS = $(wildcard *.sxpr)
DOTPS = $(subst .sxpr,.nolabel.ps,$(TRANS))
CIRCOPS = $(subst .sxpr,.nolabel.circo.ps,$(TRANS))

all: $(DOTPS) $(CIRCOPS)

%.dot %.composed: %.sxpr
	phylocomposer -d $*.dot $*.sxpr >$*.composed

%.nolabel.dot: %.dot
	cat $< | perl -pe 's/.label=.*;/;/' >$@

%.ps: %.dot
	dot -Tps $< >$@

%.circo.ps: %.dot
	circo -Tps $< >$@

%.pdf: %.ps
	ps2pdf $<

# path statistics

NFORWARD = -nf 1000000

COUNT_PATHS = grep forward-path | perl -e 'while(<>){chomp;++$$n{$$_}}for$$s(sort keys%n){print"$$n{$$s} $$s\n"}'

%.impl-stats: %.sxpr
	phylocomposer -qq $(NFORWARD) $< | $(COUNT_PATHS) >$@

%.hmmoc-stats: %.sxpr
	phylocomposer -qq $(NFORWARD) $< -ha | $(COUNT_PATHS) >$@

# don't delete secondary files
.SECONDARY:

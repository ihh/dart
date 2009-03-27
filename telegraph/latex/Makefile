TEXFILE = iterator
#TEXFILE = projects

DVIVIEWER = xdvi
PSVIEWER  = kghostview
PDFVIEWER = open

dvi: $(TEXFILE).dvi

ps:  $(TEXFILE).ps

pdf: $(TEXFILE).pdf

$(TEXFILE).dvi: $(TEXFILE).tex
	latex $(TEXFILE).tex

$(TEXFILE).ps: $(TEXFILE).dvi
	dvips -f -Ppdf -G0 $(TEXFILE).dvi >$@

$(TEXFILE).pdf: $(TEXFILE).ps
	ps2pdf $(TEXFILE).ps

bib:
	bibtex $(TEXFILE)

show: $(TEXFILE).dvi
	$(DVIVIEWER) $(TEXFILE).dvi

showps: $(TEXFILE).ps
	$(PSVIEWER) $(TEXFILE).ps

showpdf: $(TEXFILE).pdf
	$(PDFVIEWER) $(TEXFILE).pdf

clean:
	rm $(TEXFILE).dvi

# basic makefile

decon: decon.pdf

decon.pdf: decon.tex crisdefs.tex figures/*.pdf decon.bib Makefile
	pdflatex decon.tex -interaction nonstopmode && \
	bibtex decon && \
	pdflatex decon.tex -interaction nonstopmode && \
	pdflatex decon.tex -interaction nonstopmode || \
	rm decon.pdf 2> /dev/null || true

show: decon.pdf
	evince decon.pdf

install: decon.pdf
	cp -a decon.pdf ../..

clean:
	rm decon.log decon.aux decon.toc decon.vrb decon.bbl \
	decon.blg decon.snm decon.nav decon.out 2> /dev/null || true


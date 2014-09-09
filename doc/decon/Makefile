# basic makefile

airs_decon: airs_decon.pdf

airs_decon.pdf: airs_decon.tex crisdefs.tex Makefile
	pdflatex airs_decon.tex -interaction nonstopmode && \
	pdflatex airs_decon.tex -interaction nonstopmode || \
	rm airs_decon.pdf

show: airs_decon.pdf
	epdfview airs_decon.pdf

clean:
	rm airs_decon.log airs_decon.aux airs_decon.toc \
	airs_decon.snm airs_decon.nav airs_decon.out 2> /dev/null || true


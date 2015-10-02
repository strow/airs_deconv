# basic makefile

decon: decon.pdf

decon.pdf: decon.tex crisdefs.tex figures/*.pdf Makefile
	pdflatex decon.tex -interaction nonstopmode && \
	pdflatex decon.tex -interaction nonstopmode || \
	rm decon.pdf 2> /dev/null || true

show: decon.pdf
	acroread decon.pdf

clean:
	rm decon.log decon.aux decon.toc decon.vrb \
	decon.snm decon.nav decon.out 2> /dev/null || true


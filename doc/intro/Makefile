# basic makefile

intro: intro.pdf

intro.pdf: intro.tex figures/*.pdf
	pdflatex intro.tex -interaction nonstopmode || rm intro.pdf

show: intro.pdf
	acroread intro.pdf

clean:
	rm intro.log intro.aux intro.toc \
		intro.snm intro.nav intro.out 2> /dev/null || true


# basic makefile

airs2cris: airs2cris.pdf

airs2cris.pdf: airs2cris.tex crisdefs.tex figures/*.pdf
	pdflatex airs2cris.tex -interaction nonstopmode || rm airs2cris.pdf

show: airs2cris.pdf
	acroread airs2cris.pdf

clean:
	rm airs2cris.log airs2cris.aux airs2cris.toc \
	   airs2cris.snm airs2cris.nav airs2cris.out 2> /dev/null || true


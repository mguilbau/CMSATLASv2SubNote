# LaTeX Makefile
FILE=theNote

all: $(FILE).pdf view 

view:
	open $(FILE).pdf

.PHONY: clean

clean:
	\rm *.aux *.blg *.bbl *.log *.pdf

$(FILE).pdf: $(FILE).tex
	pdflatex $(FILE)
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

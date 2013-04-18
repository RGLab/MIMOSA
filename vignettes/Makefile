# Makefile to use knitr for package vignettes

# put all PDF targets here, separated by spaces
PDFS= MIMOSA.pdf

all: $(PDFS) 

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary

cleanall: clean
	rm *.pdf
%.pdf: %.Rnw
	Rscript -e "library(knitr); knit2pdf('$*.Rnw')"


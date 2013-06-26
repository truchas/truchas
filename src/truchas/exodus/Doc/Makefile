DOC= exodus-prog-guide.pdf

all: $(DOC)

clean:; rm -rf *.aux *.log *.out *.bak *.bbl *.blg
cleanest: clean; rm -rf $(DOC)

.SUFFIXES: .pdf
.tex.pdf:
	pdflatex $* && bibtex $* && pdflatex $* && pdflatex $*


TARGET=main
LATEX=pdflatex

# A4: 210 x 297 mm
# B5: 176 x 250 mm

all: ps

latexfirst:
	$(LATEX) $(TARGET)
	#CHAPTERS=`grep "\include{.*}" main.tex | sed 's/\\include{\(.*\)}/\1/'`

bibtex:
	bibtex main
# 	bibtex intro
	#bibtex Chapter1
	#bibtex Chapter2
	#bibtex Chapter3
	#bibtex Chapter4
	#bibtex Chapter5
	#bibtex Appendix1
	#bibtex Appendix2
# 	bibtex chap8
# 	bibtex chap9
# 	bibtex conclusions
latexsecond: bibtex
	$(LATEX) $(TARGET)
	makeindex $(TARGET)
	$(LATEX) $(TARGET)

latex: latexfirst latexsecond

ps: latex
	dvips -Ppdf -Pfve -t special -T 210mm,297mm -o $(TARGET).ps $(TARGET).dvi
#	dvips -Ppdf -Pfve -t a4 -o $(TARGET).ps $(TARGET).dvi

pdf: ps
	ps2pdf $(TARGET).ps

clean:
	rm -rf *.aux Appendices/*aux Chapters/*aux *.log *.bbl *.blg *.brf *.cb *.ind *.idx *.ilg  \
              *.inx *.ps *.dvi *.pdf *.toc *.out

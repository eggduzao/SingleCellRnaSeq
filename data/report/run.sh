#!/bin/bash

pdflatex report.tex
pdflatex report.tex
bibtex report
bibtex report
pdflatex report.tex
pdflatex report.tex

rm -rf report.aux
rm -rf report.bbl
rm -rf report.blg
rm -rf report.log
rm -rf report.out
rm -rf report.toc



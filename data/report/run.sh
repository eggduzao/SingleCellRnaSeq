#!/bin/bash

pdflatex report.tex -quiet > /dev/null 2>&1
pdflatex report.tex -quiet > /dev/null 2>&1
bibtex report > /dev/null 2>&1
bibtex report > /dev/null 2>&1
pdflatex report.tex -quiet > /dev/null 2>&1
pdflatex report.tex -quiet > /dev/null 2>&1

rm -rf report.aux
rm -rf report.bbl
rm -rf report.blg
rm -rf report.log
rm -rf report.out
rm -rf report.toc



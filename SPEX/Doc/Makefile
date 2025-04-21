#-------------------------------------------------------------------------------
# SPEX/Doc/Makefile
#-------------------------------------------------------------------------------

# SPEX: (c) 2019-2025, Chris Lourenco, Jinhao Chen,
# Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis.
# All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

# creates the user guide

userguide:
	pdflatex SPEX_UserGuide.tex
	bibtex SPEX_UserGuide
	pdflatex SPEX_UserGuide.tex
	pdflatex SPEX_UserGuide.tex
	pdflatex SPEX_UserGuide.tex
	- $(RM) -r *.out *.aux *.log *.bbl *.blg *.toc *.cpt

clean:
	- $(RM) -r *.out *.aux *.log *.bbl *.blg *.toc *.cpt

purge: clean

distclean: clean


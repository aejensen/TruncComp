LATEXMK = latexmk
LATEXMKFLAGS = -pdf -interaction=nonstopmode -halt-on-error
MANUSCRIPT = manuscrip.tex
MANUSCRIPT_STEM = $(basename $(MANUSCRIPT))
LATEX_TEMP_EXTENSIONS = acn acr alg aux bbl bcf blg brf fdb_latexmk fls glg glo gls idx ilg ind ist lof log lot nav out run.xml snm synctex.gz toc vrb xdv
LATEX_TEMP_FILES = $(foreach ext,$(LATEX_TEMP_EXTENSIONS),$(MANUSCRIPT_STEM).$(ext))

.PHONY: assets simulate simulate-slurm pdf all clean clean-latex-temp

assets:
	Rscript scripts/build-assets.R --target manuscript

simulate:
	Rscript ../simulation-study/scripts/run-simulation-study.R --workers=1

simulate-slurm:
	Rscript ../simulation-study/scripts/submit-simulation-study-slurm.R

pdf: assets
	$(LATEXMK) $(LATEXMKFLAGS) $(MANUSCRIPT)
	$(MAKE) clean-latex-temp

all: pdf

clean-latex-temp:
	rm -f $(LATEX_TEMP_FILES)

clean:
	$(LATEXMK) -C $(MANUSCRIPT)
	$(MAKE) clean-latex-temp
	rm -rf build

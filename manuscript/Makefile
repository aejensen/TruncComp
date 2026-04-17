LATEXMK = latexmk
LATEXMKFLAGS = -pdf -interaction=nonstopmode -halt-on-error
MANUSCRIPT = manuscript_new.tex

.PHONY: assets simulate simulate-slurm pdf all clean

assets:
	Rscript scripts/build-assets.R --target manuscript

simulate:
	Rscript ../simulation-study/scripts/run-simulation-study.R --workers=1

simulate-slurm:
	Rscript ../simulation-study/scripts/submit-simulation-study-slurm.R

pdf: assets
	$(LATEXMK) $(LATEXMKFLAGS) $(MANUSCRIPT)

all: pdf

clean:
	$(LATEXMK) -C $(MANUSCRIPT)
	rm -rf build

LATEXMK = latexmk
LATEXMKFLAGS = -pdf -interaction=nonstopmode -halt-on-error
MANUSCRIPT = manuscript.tex

.PHONY: assets simulate pdf all clean

assets:
	Rscript scripts/build-assets.R --target manuscript

simulate:
	Rscript ../simulation-study/scripts/run-simulation-study.R --workers=1

pdf: assets
	$(LATEXMK) $(LATEXMKFLAGS) $(MANUSCRIPT)

all: pdf

clean:
	$(LATEXMK) -C $(MANUSCRIPT)
	rm -rf build

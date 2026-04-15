LATEXMK = latexmk
LATEXMKFLAGS = -pdf -interaction=nonstopmode -halt-on-error
MANUSCRIPT = manuscript.tex

.PHONY: assets pdf all clean

assets:
	Rscript scripts/build-assets.R --target manuscript

pdf: assets
	$(LATEXMK) $(LATEXMKFLAGS) $(MANUSCRIPT)

all: pdf

clean:
	$(LATEXMK) -C $(MANUSCRIPT)
	rm -rf build

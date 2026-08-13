LATEXMK = latexmk
LATEXMKFLAGS = -pdf -interaction=nonstopmode -halt-on-error
MANUSCRIPT = manuscript.tex
SELF_CONTAINED_SCRIPT = scripts/build-self-contained.sh
SELF_CONTAINED_DIR = self_contained
SUBMISSION_TEX = two_component_treatment_comparisons_appendix_proofs.tex
TEX_SOURCES = $(MANUSCRIPT)
TEX_STEMS = $(basename $(TEX_SOURCES))
LATEX_TEMP_EXTENSIONS = acn acr alg aux bbl bcf blg brf dvi fdb_latexmk fls glg glo gls idx ilg ind ist loa lof log lol lom lot maf mtc mtc0 nav nlo nls out ps run.xml snm spl synctex.gz thm tmp toc vrb xdv
LATEX_TEMP_FILES = $(foreach stem,$(TEX_STEMS),$(foreach ext,$(LATEX_TEMP_EXTENSIONS),$(stem).$(ext)))
LATEX_GENERATED_TEX = $(wildcard *-pratend*.tex)

.NOTPARALLEL:
.PHONY: assets simulate simulate-slurm pdf self-contained submission-pdf all clean clean-latex-temp

assets:
	Rscript scripts/build-assets.R --target manuscript

simulate:
	Rscript ../simulation-study/scripts/run-simulation-study.R --workers=1

simulate-slurm:
	Rscript ../simulation-study/scripts/submit-simulation-study-slurm.R

pdf: assets
	@trap '$(MAKE) --no-print-directory clean-latex-temp' 0; \
	$(LATEXMK) $(LATEXMKFLAGS) $(MANUSCRIPT)

self-contained:
	./$(SELF_CONTAINED_SCRIPT)

submission-pdf:
	@cd "$(SELF_CONTAINED_DIR)" && \
	trap '$(LATEXMK) -c "$(SUBMISSION_TEX)" >/dev/null 2>&1 || true; $(RM) $(basename $(SUBMISSION_TEX))-pratend*.tex' 0; \
	$(LATEXMK) $(LATEXMKFLAGS) "$(SUBMISSION_TEX)"

all: pdf

clean-latex-temp:
	rm -f $(LATEX_TEMP_FILES) $(LATEX_GENERATED_TEX)

clean: clean-latex-temp

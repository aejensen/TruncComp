# Rewrite Abstract and Introduction Prompt

Use this prompt when you want an agent to rewrite only the abstract and
Introduction of `manuscript/manuscrip.tex`, with the randomized-trial
single-p-value motivation made central.

```text
Write a new abstract and Introduction for the TruncComp manuscript.

Use manuscript/manuscrip.tex as the manuscript file. If the task says manuscript.tex, treat that as referring to manuscript/manuscrip.tex.

Before writing, first read:
Kurland BF, Johnson LL, Egleston BL, Diehr PH. "Longitudinal Data with Follow-up Truncated by Death: Match the Analysis Method to Research Aims." Statistical Science. 2009;24(2):211. DOI: 10.1214/09-STS293.
Accessible source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2812934/

Use Kurland et al. to position this manuscript. Their taxonomy distinguishes unconditional models, fully conditional models, principal-stratification/causal models, partly conditional survivor models, and joint models of survival and longitudinal response. Explain where this manuscript fits: it is not a survivor-only or principal-stratum analysis; it targets an observed-data atom-plus-non-atom endpoint at a fixed follow-up time. It is closest in spirit to an observed-data joint/composite-response aim, but sharper than a one-dimensional coded composite because it preserves the atom/non-atom decomposition while still producing one joint test.

Read the full manuscript first, especially the Method section, Bayesian two-group comparison subsection, Simulation study, IST-3 application, and Discussion.

Follow AGENTS.md. Do not rerun, refresh, replace, or invalidate manuscript/application-data/ist3-bayes-cache.rds. Do not regenerate generated tables, figures, diagnostics, Bayesian outputs, application assets, or simulation assets.

Goal:
Rewrite the abstract and Introduction so they are publication-ready, concise, and aligned with the manuscript's actual contribution, with special emphasis on randomized controlled trials where a single prespecified primary p-value is often desirable.

The new abstract should:
1. State the RCT problem: continuous patient-centered outcomes can be undefined after death or another substantive atom, yet trials often require one primary analysis and one p-value.
2. Explain why coded-scale or rank-based one-dimensional analyses provide a single comparison but can conflate atom-probability and non-atom outcome effects.
3. State the main contribution: parametric and semi-parametric likelihood-ratio tests that give one prespecified p-value for the joint observed-data null of no treatment effect on the atom probability and non-atom mean component.
4. Mention that the framework also reports interpretable component summaries and confidence regions, so the single p-value does not come at the cost of hiding the two components.
5. Briefly mention the Bayesian nonparametric two-group counterpart using flexible mixture priors.
6. Summarize simulations and the IST-3 RCT application without causal or clinical overclaiming.

The new Introduction should:
1. Motivate truncation by death and atom-plus-continuous endpoints in RCTs, where regulators, protocols, and trial reports often favor a single prespecified primary p-value.
2. Use Kurland et al. to explain that analysis methods should match research aims, then place this manuscript in their taxonomy as an observed-data atom-plus-non-atom endpoint analysis at fixed follow-up.
3. Distinguish this aim from survivor-only analyses, partly conditional survivor summaries, and principal-stratum causal estimands.
4. Explain that one-dimensional composite analyses answer a valid coded/ranked endpoint question, but can blur whether treatment acts through survival/atom probability, the non-atom outcome, or both.
5. Introduce the proposed two-component joint null as a way to retain one primary p-value while preserving interpretable component-level summaries.
6. Preview the frequentist parametric and semi-parametric LRTs.
7. Preview the Bayesian nonparametric two-group comparison as a complementary model-based analysis of arm-specific observed-data laws, not as a replacement for the frequentist single-p-value test.
8. Preview the simulations, IST-3 application, and discussion.

Style constraints:
- Preserve LaTeX style and existing citation commands.
- Add the Kurland et al. citation to bibliography only if it is not already present.
- Do not overclaim robustness, efficiency, causal interpretation, regulatory acceptance, or clinical benefit.
- Distinguish frequentist semi-parametric empirical likelihood from Bayesian nonparametric mixture modeling.
- Keep the Introduction close to the current length unless modest expansion improves clarity.
- Keep the RCT/single-p-value motivation practical and central, but do not make the paper sound like a regulatory guidance document.

Output:
1. Briefly diagnose the current abstract/introduction's main weaknesses.
2. Provide complete replacement LaTeX blocks for the abstract environment and \section{Introduction}.
3. If editing directly, only edit the abstract and Introduction in manuscript/manuscrip.tex, plus bibliography.bib only if the Kurland citation is missing.
4. Verify with:
   cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error manuscrip.tex
   Do not use make pdf.
5. Report whether compilation succeeded and confirm the IST-3 Bayesian cache and generated assets were not touched.
```

## Notes

- The RCT single-p-value motivation should be central in both the abstract and
  Introduction.
- The Kurland taxonomy should position the research aim, not turn the manuscript
  into a longitudinal trajectory paper.

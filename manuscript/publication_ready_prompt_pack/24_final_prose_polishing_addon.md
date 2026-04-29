# Add-On Prompt for Final Prose Polishing

```text
Use this add-on for phases 13-17 of the TruncComp manuscript.

Preserve technical precision while improving readability. For each paragraph,
identify its function:
- motivation;
- definition;
- assumption;
- construction;
- formal result;
- proof explanation;
- implementation;
- numerical evidence;
- Bayesian model or posterior summary;
- limitation;
- transition.

Split paragraphs that perform incompatible functions. Merge or shorten
redundant paragraphs. Do not add broad claims, clinical efficacy claims,
principal-stratum causal claims, or robustness claims not supported by the
manuscript.

Every section should make clear:
1. Why it is needed.
2. What object, result, evidence, or limitation it introduces.
3. How it advances the two-component observed-data argument.
4. How it prepares the next section.

For final verification, use direct latexmk on manuscrip.tex if verification is
authorized. Do not use make pdf for text-only edits.
```

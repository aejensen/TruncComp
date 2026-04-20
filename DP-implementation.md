# DP Implementation Memo

This document captures implementation ideas and preliminary conclusions from our discussion.

## Entry Template

Each logged entry includes:

- **Name**
- **Idea Description**
- **Commentary / Preliminary Conclusion**

---

## Log

### Bayesian Two-Part Formulation for `TruncComp2`

- **Name**: Bayesian two-part model for the observed indicator and continuous observed outcome
- **Idea Description**: Formulate the outcome in two pieces: an indicator `A = 1(Y != atom)` for whether the continuous outcome is observed, and the continuous value `Y` conditional on `A = 1`. Model the binary part with Bayesian logistic regression, for example `A_i | R_i, L_i ~ Bernoulli(pi_i)` with `logit(pi_i) = alpha_0 + alpha_R R_i + f_alpha(L_i)`. Model the continuous part on the observed rows only, for example through a Bayesian regression for `Y_i | A_i = 1, R_i, L_i`. A derived combined estimand can be defined through the extended mean `E[Y_ext | R, L] = atom * (1 - pi_R(L)) + pi_R(L) * m_R(L)`, where `m_R(L)` is the conditional mean among observed outcomes. Posterior draws can then be used to summarize treatment effects on the observation probability, the observed continuous outcome, and the combined extended outcome.
- **Commentary / Preliminary Conclusion**: The main shift from the current `SPLRT` framing is to replace likelihood-ratio testing with posterior inference. A practical first implementation would use Bayesian logistic regression for the observation model together with a Bayesian normal or Student-`t` regression for the observed continuous outcome, and then derive posterior draws for the combined estimand. If a more genuinely semi-parametric Bayesian version is desired, the continuous component could later be extended to a Dirichlet process mixture of normals, with BART or Gaussian-process style models as more flexible but less low-dimensional alternatives. For an initial package implementation, posterior means, credible intervals, and posterior probabilities such as `Pr(Delta > 0 | data)` seem more robust and interpretable than Bayes factors.

### Bayesian Single-Number Summary for Group Comparison

- **Name**: How to summarize a Bayesian group comparison by a single number
- **Idea Description**: Let `Delta` denote the target group-comparison estimand, ideally the treatment-group difference in the combined extended outcome. In a Bayesian analysis there is no exact analogue of the frequentist p-value that is universally preferred. The most direct single-number summary is the posterior probability of benefit, for example `Pr(Delta > 0 | data)` when larger values are better. If the scientific question is about a clinically relevant improvement rather than merely a positive difference, a more decision-oriented summary is `Pr(Delta > delta_clin | data)` for a pre-specified threshold `delta_clin`. Another possibility is a Bayes factor or posterior odds comparing `H0: Delta = 0` to an alternative model, but that requires an explicit prior on both the null and the alternative and behaves differently from standard posterior summaries.
- **Commentary / Preliminary Conclusion**: The cleanest default for the memo seems to be that the Bayesian alternative to a single p-value is not a Bayes factor but a posterior probability attached to the estimand of interest. If the goal is a directional summary, report `Pr(Delta > 0 | data)`. If the goal is decision support, report `Pr(Delta > delta_clin | data)` for a meaningful threshold. This keeps the summary tied directly to the treatment effect that users care about and avoids the prior sensitivity of Bayes factors for point-null testing. Bayes factors may still be useful as an optional secondary output when the scientific aim is explicit model comparison, but they should probably not be the default headline summary.

### Different Ways to Use a Dirichlet Process

- **Name**: Possibilities for introducing a Dirichlet process in the Bayesian formulation
- **Idea Description**: There are several distinct places where a Dirichlet process can be used in the two-part model, and they are not equivalent. One option is a Dirichlet process mixture for the observed continuous outcome model, for example by writing `Y_i | A_i = 1, R_i, L_i` as a regression mean plus an error term whose distribution is modeled by a DP mixture of normals. This is the most direct semiparametric extension because it preserves a simple treatment-effect parameter in the mean model while allowing the observed outcome distribution to be flexible. A second option is to place a DP prior directly on the conditional distribution of the observed outcome within each treatment arm, possibly separately by arm. That gives great flexibility for density estimation but makes it less natural to define a single regression-style treatment parameter when covariates are present. A third option is to use a dependent Dirichlet process, where the random mixing distribution varies with treatment or covariates. That is more expressive and can capture treatment-specific or covariate-specific shape changes in the outcome distribution, but it is also more complex computationally and conceptually. A fourth possibility is to use a DP only as a random-effects or latent-clustering device, allowing subjects to be grouped into latent subclasses with different outcome distributions; this can be useful if one expects hidden heterogeneity beyond observed covariates.
- **Commentary / Preliminary Conclusion**: The most attractive first use of a Dirichlet process is in the continuous component only, through a DP mixture on the residual or observed-outcome distribution conditional on `A = 1`. That keeps the binary observation model simple, preserves interpretable treatment-effect summaries, and is closest in spirit to a semiparametric replacement for the current frequentist continuous-component model. Using separate arm-specific DPs is appealing for flexibility but risks making group comparisons harder to summarize by one parameter. Dependent DPs are scientifically interesting if treatment is expected to change the shape of the observed outcome distribution and not just its mean, but they may be too ambitious for an initial implementation in `TruncComp2`. A reasonable staged plan would be to start with parametric Bayesian logistic plus a DP-mixture continuous model, and only later consider dependent-DP formulations if there is a strong need for distributional treatment effects beyond the mean.

### Separate DP Models by Treatment Group

- **Name**: Pros and cons of fitting separate Dirichlet process models in each treatment arm
- **Idea Description**: One especially flexible strategy is to model the observed continuous outcome distribution separately within each treatment group, for example by assigning independent or partially linked Dirichlet process priors to the two arm-specific distributions for `Y | A = 1, R = 0` and `Y | A = 1, R = 1`. This avoids imposing a shared parametric form or a simple location-shift structure across treatment arms. Each group can have its own skewness, tail behavior, multimodality, and local features, and treatment effects can then be defined afterward as posterior contrasts of functionals such as mean differences, median differences, quantile differences, probabilities of superiority, or contrasts in the combined extended-outcome mean.
- **Commentary / Preliminary Conclusion**: The major strength of separate arm-specific DP models is that they let treatment affect the full observed-outcome distribution in whatever way the data support, which is scientifically appealing when one does not want to prejudge the form of the effect. They are also robust to misspecification and naturally support rich posterior contrasts beyond the mean. The main cost is that they do not produce a single built-in treatment coefficient; the treatment effect must instead be defined through a chosen functional of the two posterior distributions. They also remove borrowing across treatment groups unless a higher-level hierarchical structure is added, which may hurt efficiency when one arm has few observed outcomes. Covariate adjustment is another challenge, because once baseline covariates are included the model usually needs to move beyond two separate unconditional DPs toward conditional or dependent-DP constructions. Overall, separate DP models by treatment arm are a strong choice if the priority is maximal flexibility in the outcome distribution, but they are less attractive if the priority is a simple regression-style estimand, straightforward adjustment, or a default one-number summary. A useful compromise may be to consider hierarchical arm-specific DP priors later, so that the groups can differ flexibly while still sharing information when appropriate.

### Two Ways to Model the Response Within Each Treatment Arm

- **Name**: Compare explicit Bernoulli-plus-DP mixtures with augmented-space DP mixtures
- **Idea Description**: With separate models in each treatment arm and no covariates, two natural formulations arise. In the first, each arm `r` is modeled explicitly as a two-part mixture
  `Y | R = r ~ p_r delta_atom + (1 - p_r) F_r^c`,
  where `p_r` is the probability of the atom and `F_r^c` is a Dirichlet process mixture on the continuous support away from the atom. In the second, each arm is given one unified Dirichlet process mixture on an augmented support that includes both a point-mass component at the atom and continuous components elsewhere. A convenient way to write that is
  `Y | R = r, G_r ~ int K(y | theta) G_r(d theta)`,
  where the kernel family contains one special kernel equal to `delta_atom` and the rest are continuous kernels, and `G_r` has a DP prior with base measure that places positive mass on both the atom component and the continuous part.
- **Commentary / Preliminary Conclusion**: The explicit Bernoulli-plus-DP formulation is easier to interpret because it keeps the atom probability `p_r` and the continuous response distribution `F_r^c` as separate primary objects. That matches the scientific decomposition already used in the frequentist formulation, makes prior specification simpler because one can put a direct prior on `p_r`, and gives straightforward posterior summaries such as differences in atom probability, differences in the mean among non-atom outcomes, and differences in the combined extended mean. It is also computationally appealing because the atom part and the continuous part can be updated in a modular way. The augmented-space DP mixture is more elegant as a single unified model for the full response distribution and is attractive if the main goal is to estimate each treatment-arm distribution nonparametrically without committing in advance to a two-part interpretation. However, it gives less direct control over the prior on the atom probability, since the mass at the atom is induced by the DP prior rather than parameterized separately, and it can make interpretation slightly less transparent because the atom and continuous pieces are no longer foregrounded as separate parameters. In practice, the two approaches may be closer than they first appear: if the augmented model uses an exact point-mass kernel at the atom and continuous kernels elsewhere, then it effectively induces an atom probability together with a continuous conditional distribution, so it can often be viewed as a more tightly coupled reparameterization of the explicit two-part model. For this application, the explicit Bernoulli-plus-DP construction seems preferable if the goal is interpretability and continuity with the current `TruncComp` logic, whereas the augmented-space formulation is preferable if the goal is to present the analysis as a single fully nonparametric model for the entire arm-specific response distribution.

### Measure-Theoretic Formulation and Equivalence of the Two Response Models

- **Name**: Formal write-up of the explicit two-part model and the augmented-space DP model
- **Idea Description**: Let `(S, \mathcal{B})` be the measurable outcome space, let `a \in S` denote the special atom, and let `C = S \setminus \{a\}`. Let `(Theta_c, \mathcal{T}_c)` be a measurable parameter space for continuous kernels, and let `K_c(· | theta)` be a Markov kernel on `(S, \mathcal{B})` such that `K_c({a} | theta) = 0` for every `theta \in Theta_c`. This encodes that the continuous part places no mass on the atom. For a fixed treatment arm `r`, the explicit two-part model can then be written as

```math
p_r \in [0,1], \qquad G_r^c \in \mathcal{P}(\Theta_c),
```

```math
F_r(B) = p_r \delta_a(B) + (1 - p_r) \int_{Theta_c} K_c(B \mid \theta)\, G_r^c(d\theta),
\qquad B \in \mathcal{B},
```

with observations satisfying `Y_{ri} | F_r iid F_r`. A natural Bayesian version is

```math
p_r \sim \Pi_r, \qquad G_r^c \sim DP(\gamma_r H_r),
```

where `H_r` is a probability measure on `(Theta_c, \mathcal{T}_c)`. The augmented-space formulation starts from the enlarged parameter space `Theta = \{*\} \cup Theta_c`, where `*` is a formal symbol representing the atom kernel. Define the combined kernel

```math
K(B \mid *) = \delta_a(B), \qquad K(B \mid \theta) = K_c(B \mid \theta), \ \theta \in Theta_c.
```

Let

```math
\widetilde{H}_r = w_r \delta_* + (1 - w_r) H_r,
\qquad 0 < w_r < 1,
```

and place a Dirichlet process prior on the full mixing measure:

```math
G_r \sim DP(\alpha_r \widetilde{H}_r),
\qquad
F_r(B) = \int_{Theta} K(B \mid \theta)\, G_r(d\theta).
```

The equivalence is exact when the explicit model uses the matching prior

```math
p_r \sim Beta(\alpha_r w_r, \alpha_r (1 - w_r)),
\qquad
G_r^c \sim DP(\alpha_r (1 - w_r) H_r),
```

with `p_r` independent of `G_r^c`. In that case define a random probability measure on `Theta` by

```math
G_r = p_r \delta_* + (1 - p_r) G_r^c.
```

Then for every measurable partition `B_1, ..., B_k` of `Theta_c`,

```math
\bigl(G_r({*}), G_r(B_1), ..., G_r(B_k)\bigr)
= \bigl(p_r, (1-p_r)G_r^c(B_1), ..., (1-p_r)G_r^c(B_k)\bigr)
```

has Dirichlet distribution

```math
Dir\bigl(\alpha_r w_r,\ \alpha_r (1-w_r) H_r(B_1),\ ..., \ \alpha_r (1-w_r) H_r(B_k)\bigr),
```

which is precisely the finite-dimensional law of `DP(\alpha_r \widetilde{H}_r)`. Hence `G_r \sim DP(\alpha_r \widetilde{H}_r)`, and the induced response law becomes

```math
F_r(B)
= \int_{Theta} K(B \mid \theta)\, G_r(d\theta)
= p_r \delta_a(B) + (1-p_r)\int_{Theta_c} K_c(B \mid \theta)\, G_r^c(d\theta),
```

which is exactly the explicit two-part model. Conversely, if `G_r \sim DP(\alpha_r \widetilde{H}_r)` in the augmented model, define

```math
p_r = G_r({*}),
\qquad
G_r^c(A) = \frac{G_r(A)}{1-p_r}, \quad A \subseteq Theta_c,
```

on the event `p_r < 1`. Then the restriction property of the Dirichlet process yields

```math
p_r \sim Beta(\alpha_r w_r, \alpha_r(1-w_r)),
\qquad
G_r^c \sim DP(\alpha_r(1-w_r) H_r),
```

with `p_r` independent of `G_r^c`, and the same decomposition of `F_r` follows.
- **Commentary / Preliminary Conclusion**: In this sense the two models are not merely similar but equivalent parameterizations of the same prior on the arm-specific response distribution, provided that the augmented model contains an exact point-mass kernel at the atom, the remaining kernels assign zero probability to the atom, and the explicit two-part prior on `p_r` is chosen to match the Beta law induced by the augmented-space DP. The practical difference is therefore not the set of response distributions being modeled, but the parameterization and the transparency of the decomposition. The equivalence breaks if one wants an arbitrary prior on `p_r` rather than the Beta law implied by the DP, or if the continuous kernels are allowed to place nonzero mass on the atom. Under those changes the augmented-space model becomes strictly more general, but also less cleanly interpretable as a separated atom-plus-continuous response model.

### Modeling and Interpretation Differences Between the Two Approaches

- **Name**: Modeling and interpretation contrast between the explicit two-part model and the augmented-space DP model
- **Idea Description**: Even when the two formulations are formally equivalent under matched priors, they emphasize different scientific stories. The explicit two-part model starts from the view that the atom and the continuous responses are qualitatively different phenomena and should therefore be modeled separately. In that formulation, the atom probability and the continuous non-atom distribution are primitive parameters, so the model naturally supports statements such as “treatment changes the probability of the atom” and “among non-atom outcomes, treatment changes the response distribution.” The augmented-space DP model starts instead from the view that each treatment arm has one unknown response distribution on the full outcome space, and the atom is simply one part of that distribution. In that formulation, the primary object is the full arm-specific random probability measure, and the decomposition into atom mass and conditional continuous distribution is secondary or induced after the fact.
- **Commentary / Preliminary Conclusion**: From a modeling perspective, the explicit two-part approach is better aligned with the original `TruncComp` logic because it builds the special role of the atom directly into the model structure. That makes prior elicitation easier whenever investigators have different beliefs about the atom frequency and the shape of the observed continuous outcomes, and it makes derived estimands easier to explain because the decomposition is already built into the generative story. The augmented-space formulation is conceptually cleaner if the goal is to say only that each treatment arm has an unknown distribution and to avoid privileging the atom structurally beyond giving it a point-mass kernel. This can feel more “fully nonparametric,” but it may also hide the scientific distinction between the atom event and the continuous response unless that distinction is reintroduced in posterior summaries. So the main difference is not mathematical flexibility so much as which interpretation is primary: the explicit model treats the problem as one of decomposing two substantively different mechanisms, while the augmented model treats it as one of nonparametric estimation of a single arm-specific distribution. For this project, the explicit two-part framing still seems preferable if interpretability and continuity with the current methodology are central, whereas the augmented-space formulation is most attractive if one wants the arm-specific distribution itself to be the main inferential object.

### Implementation and Fitting Differences Between the Two Approaches

- **Name**: Implementation and computation contrast between the explicit two-part model and the augmented-space DP model
- **Idea Description**: From a fitting perspective, the explicit two-part model naturally decomposes into two subproblems within each treatment arm: estimate the atom probability from all outcomes, and fit the DP mixture only to the non-atom observations. With a conjugate prior, the atom probability can often be updated in closed form through a Beta-Bernoulli step, while the continuous part can be fit using standard DP-mixture machinery on the reduced dataset consisting only of observed continuous responses. By contrast, the augmented-space DP model treats all observations jointly inside one mixture model, with the atom represented by a distinguished kernel or component in the mixture. This can be elegant if one wants a single sampler and a single representation of the arm-specific distribution, but it typically requires the fitting algorithm to handle a special exact point-mass component alongside continuous kernels.
- **Commentary / Preliminary Conclusion**: In implementation terms, the explicit two-part formulation is usually easier to code, debug, and validate. It is modular, because the Bernoulli part and the continuous DP-mixture part can be developed and tested separately, and it is often computationally cheaper because the continuous sampler only runs on non-atom observations. Posterior summaries are also direct: the atom probability is an explicit parameter, and the continuous posterior is already conditioned on being away from the atom. The augmented-space formulation has the advantage of giving one unified inference engine for the entire response distribution, which can be attractive if the software goal is to expose a single nonparametric object rather than two linked components. However, it is often harder to implement with off-the-shelf DP-mixture code, since many standard routines assume only continuous kernels and do not naturally accommodate an exact Dirac component. It can also create more subtle mixing and identifiability issues if the continuous kernels can place appreciable mass near the atom, because the sampler then has to distinguish “true atom mass” from “continuous clusters concentrated close to the atom.” Even when the two models are formally equivalent under matched priors, the explicit parameterization is likely to behave better computationally and to be simpler to maintain in a package. The augmented-space formulation becomes more attractive when one is willing to invest in custom samplers or when the conceptual benefit of a single unified random distribution is worth the extra implementation complexity.

### Recommended DP Mixture Structure for the First Implementation

- **Name**: Recommended structure for the no-covariate arm-specific Bayesian DP model
- **Idea Description**: Based on the discussion so far, the recommended first construction is to fit separate arm-specific models and to use the explicit two-part parameterization within each arm. For each treatment arm `r in {0,1}`, write

```math
Y \mid R = r \sim p_r \delta_a + (1-p_r) F_r^c,
```

where `a` is the atom, `p_r` is the atom probability in arm `r`, and `F_r^c` is the continuous non-atom response distribution. Model `p_r` directly, for example with a Beta prior,

```math
p_r \sim Beta(a_p, b_p),
```

and model `F_r^c` as a Dirichlet process mixture on the support `C = S \setminus \{a\}`:

```math
G_r^c \sim DP(\alpha_r H_r),
\qquad
F_r^c(B) = \int_{\Theta_c} K_c(B \mid \theta)\, G_r^c(d\theta),
\quad B \subseteq C.
```

The kernel family `K_c` should be chosen to match the support of the observed continuous outcome, for example Gaussian kernels on `\mathbb{R}`, log-normal or gamma-style kernels on `(0,\infty)`, or transformed kernels on a closed interval. The same kernel family should be used in both treatment arms so that posterior contrasts remain comparable, although the mixing distributions `G_0^c` and `G_1^c` are allowed to differ freely. Treatment effects should then be defined as posterior contrasts of functionals of the two fitted arm-specific distributions, especially the atom-probability difference, the mean difference among non-atom outcomes, and the difference in the combined extended mean

```math
\Delta = \left[a p_1 + (1-p_1) \mu_1^c\right] - \left[a p_0 + (1-p_0) \mu_0^c\right],
\qquad
\mu_r^c = \int_C y\, F_r^c(dy).
```

- **Commentary / Preliminary Conclusion**: This seems like the right recommendation for the first implementation because it preserves the substantive atom-versus-continuous decomposition, matches the existing `TruncComp` intuition, and is the cleanest computationally. It also keeps open the option of using the augmented-space formulation later, since the two formulations can be matched exactly under suitable priors, but it does not force the package to start with the more difficult unified parameterization. The main practical recommendation is therefore: separate models by treatment arm, explicit Bernoulli-plus-DP structure within arm, and support-appropriate continuous kernels that assign zero mass to the atom. A sensible default is to keep the two arm-specific priors independent in the first version for simplicity, while leaving hierarchical borrowing across arms as a later enhancement if small-sample efficiency becomes a concern. In short, the recommended architecture is not “one DP over everything,” but rather “one explicit atom model plus one continuous DP mixture per arm,” with posterior group comparison performed through functionals of the fitted arm-specific response distributions.

### Is a DP Mixture the Right Model in This Restricted Setting?

- **Name**: Recommendation on whether a Dirichlet process mixture is the most appropriate model here
- **Idea Description**: In the restricted setting currently under discussion, there are no covariates and the response is represented as a special atom together with a continuous non-atom component, modeled separately within each treatment arm. The question is whether the continuous non-atom distribution should be modeled with a Dirichlet process mixture or whether some simpler or different Bayesian construction would be preferable.
- **Commentary / Preliminary Conclusion**: The current conclusion is yes: a Dirichlet process mixture appears to be the most appropriate default model for the continuous non-atom part in this restricted setting, provided it is embedded inside the explicit two-part construction rather than used as a single undifferentiated model for the whole outcome. This recommendation follows from the combination of flexibility, robustness to distributional misspecification, compatibility with separate arm-specific modeling, and continuity with the overall `TruncComp` decomposition into atom versus observed continuous outcome. It is not meant as a universal claim that a DP mixture is always optimal. If the sample size is very small, the DP mixture may be more elaborate than the data can support, and if the inferential goal is only a mean contrast then a simpler robust parametric model could be sufficient. The kernel and base measure must also be chosen to respect the support of the observed outcome. Still, for the present restricted case, the recommended default remains: separate arm-specific models, explicit modeling of the atom probability, and a DP mixture for the continuous non-atom response distribution.

### Adding Covariates to the Bayesian Two-Part Model

- **Name**: Possibilities for incorporating covariates in the Bayesian model
- **Idea Description**: Once covariates `L` are introduced, the natural extension is still to preserve the explicit two-part structure, but now allow both the atom probability and the continuous non-atom distribution to depend on `L`. At the highest level this means specifying, for each treatment arm `r`, a conditional atom probability `p_r(L) = P(Y = a \mid R = r, L)` and a conditional continuous distribution `F_r^c(\cdot \mid L)` on the non-atom support. There are several ways to do this depending on the type and dimension of the covariates. For low-dimensional categorical covariates, one option is stratification: fit separate two-part models within each covariate cell, possibly with hierarchical priors to borrow strength across sparse cells. This is simple and interpretable, but becomes unstable when there are many categories or many empty combinations. For continuous covariates, a more practical option is to model the atom part with a regression such as logistic or probit regression and to model the continuous part through a regression-plus-DP construction, for example `Y = m_r(L) + epsilon` among non-atom observations with the residual distribution modeled by a DP mixture. This allows covariates to affect the mean or location through `m_r(L)` while the DP handles non-normality, multimodality, and heavy tails. For a mix of categorical and continuous covariates, one can combine the two ideas by using indicator terms, smooth terms, or additive functions inside the atom model and the continuous mean model, with a DP mixture for the residual distribution. A more ambitious alternative is to let the full mixing distribution depend on covariates, that is, to use a dependent Dirichlet process or predictor-dependent mixture model so that the entire conditional outcome distribution may change with `L` and not only its mean. This includes constructions such as ANOVA-DDP for categorical covariates, kernel stick-breaking, probit stick-breaking, or other covariate-dependent mixture-weight schemes for continuous or mixed covariates.
- **Commentary / Preliminary Conclusion**: The main practical distinction is between using covariates to shift a relatively simple regression structure and using covariates to change the whole mixture distribution. For categorical covariates with only a few levels, stratified or hierarchically pooled arm-specific models are attractive because they preserve interpretability and stay close to the no-covariate formulation. For continuous covariates, and especially for a mix of continuous and categorical covariates, the most sensible first extension seems to be regression in both parts of the model: a logistic or probit model for the atom probability and a regression-plus-DP-mixture model for the continuous non-atom outcome. This is much easier to fit than a fully dependent DP and still gives substantial robustness. Fully covariate-dependent DPs are scientifically appealing if one wants covariates to affect the full conditional distribution in a highly flexible way, but they introduce much more computational and modeling complexity and may be hard to stabilize in a first package implementation. A good staged recommendation is therefore: use explicit two-part modeling throughout; for categorical covariates, consider stratification with hierarchical pooling if dimensions are small; for continuous or mixed covariates, use regression on the atom probability and on the continuous component together with a DP mixture for residuals; and reserve fully dependent-DP constructions for later versions when there is a clear need for richer distributional covariate effects. In all cases, if the target is a marginal treatment contrast rather than a conditional one, posterior summaries should be standardized over the empirical or a chosen reference distribution of the covariates.

### Additive Regression Plus DP-Mixture Residuals

- **Name**: Explanation of additive regression with Dirichlet process mixture residuals
- **Idea Description**: The basic idea is to separate systematic covariate effects from unexplained distributional shape. In the continuous non-atom part of the model, one writes the response among non-atom observations as

```math
Y_i = m_r(L_i) + \varepsilon_i, \qquad \text{for } Y_i \neq a \text{ and } R_i = r,
```

where `m_r(L)` is a structured regression function and `\varepsilon_i` is a residual term whose distribution is modeled flexibly by a Dirichlet process mixture. An additive regression specification means that the regression function is built from a sum of lower-dimensional covariate contributions rather than a fully unrestricted surface, for example

```math
m_r(L) = \beta_{0r} + \sum_{j=1}^p f_{jr}(L_j),
```

where each `f_{jr}` may be a linear term, a spline, a step function for a categorical covariate, or another low-dimensional smooth effect. The residual distribution is then modeled as

```math
\varepsilon_i \mid G_r \sim \int K(\cdot \mid \theta)\, G_r(d\theta),
\qquad
G_r \sim DP(\alpha_r H_r),
```

so that the errors may be skewed, heavy-tailed, or multimodal without having to specify a single parametric family. Equivalently, one may write

```math
Y_i \mid L_i, R_i=r, G_r \sim \int K(\cdot \mid m_r(L_i), \theta)\, G_r(d\theta),
```

where the kernel is centered at the additive mean function and the DP mixture acts on the remaining unexplained variation. In the atom component, the analogous additive idea would be

```math
\text{logit}\, p_r(L) = \alpha_{0r} + \sum_{j=1}^p g_{jr}(L_j),
```

so the overall two-part model remains additive in covariates while allowing flexible residual behavior in the continuous part.
- **Commentary / Preliminary Conclusion**: This construction is attractive because it balances interpretability and flexibility. The additive regression part captures the main systematic covariate effects in a way that remains understandable and estimable even when covariates are mixed or moderately high-dimensional, while the DP mixture is reserved for what is left over after those mean effects are removed. That avoids asking the DP to learn the entire covariate-response surface nonparametrically, which would be much harder computationally and statistically. It also gives a clear interpretation to covariate effects: the functions `f_{jr}` describe how the expected non-atom response changes with each covariate, while the DP mixture accounts for deviations from simple Gaussian residual assumptions. In practice this seems like the most useful middle ground for a first covariate-adjusted Bayesian implementation. It is richer than ordinary additive regression because it allows nonstandard residual distributions, but much simpler than a fully predictor-dependent DP. One caveat is that if covariates affect not just the mean but also the residual shape or variance in strong ways, a pure additive mean plus common residual-mixture model may still be too restrictive, in which case one might later allow covariate-dependent residual scales or treatment- and covariate-dependent mixture weights.

### Special Case: Only a Categorical Covariate

- **Name**: Recommended handling when covariates are purely categorical
- **Idea Description**: If the covariate structure is purely categorical, the problem becomes much simpler because one can condition directly on covariate levels rather than building a functional regression model. In that case, for each treatment arm `r` and covariate level `l`, a natural model is

```math
Y \mid R=r, L=l \sim p_{rl}\delta_a + (1-p_{rl})F^c_{rl},
```

where `p_{rl}` is the atom probability in the `(r,l)` cell and `F^c_{rl}` is the corresponding continuous non-atom response distribution, for example modeled by a DP mixture. This yields a fully stratified formulation in which each treatment-covariate cell has its own two-part model. A more stable extension is to keep the same cell-specific decomposition but add hierarchical priors that partially pool the `p_{rl}` values and/or the continuous-distribution parameters across levels of `l`. An alternative is an ANOVA-style regression parameterization for the atom part and for the continuous mean, but for purely categorical covariates this may be less direct than working with cell-specific quantities.
- **Commentary / Preliminary Conclusion**: For categorical covariates only, the most attractive first strategy seems to be direct stratification, because it preserves the interpretation of the no-covariate model and makes the effect of the covariate easy to explain. If the number of categories is small and the cells are well populated, fully separate cell-specific two-part models are likely to be the simplest and clearest approach. If some cells are sparse, hierarchical pooling across categories is preferable so that the model can borrow strength while retaining cell-specific interpretation. In this setting, additive regression is probably not the first tool to reach for unless there is a strong reason to impose a lower-dimensional ANOVA-type structure. So the practical recommendation for purely categorical covariates is: start with stratified explicit two-part models, add hierarchical pooling if needed, and only move to more formal regression parameterizations if the number of categories or interactions makes direct stratification unwieldy.

### Conditional Effects With Categorical Covariates in the DP Model

- **Name**: How to recover a frequentist-style adjusted conditional effect when `L` is categorical
- **Idea Description**: In the frequentist adjusted model, the reported treatment effect is conditional on `L` and is common across covariate strata under a no-interaction assumption. The same distinction appears in the Bayesian DP setting. For the continuous non-atom part, one can certainly stratify by treatment arm and allow the mean to depend on the categorical covariate:

```math
Y_i \mid (Y_i \neq a, R_i=r, L_i=l, G_r)
\sim \int K(\cdot \mid m_r(l), \theta)\, G_r(d\theta),
```

where `m_r(l)` is a level-specific mean function for arm `r`. In this formulation the conditional treatment effect is naturally

```math
\mu(l) = m_1(l) - m_0(l),
```

so the treatment effect is stratum-specific. If the goal is instead to mimic the adjusted frequentist estimand, then one must impose a common-shift restriction across levels of `L`, for example

```math
m_r(l) = \beta_0 + \mu r + \beta_l,
```

or equivalently

```math
m_1(l) - m_0(l) = \mu \qquad \text{for all } l.
```

The same logic applies to the atom part. A fully stratified model with cell-specific atom probabilities `p_{rl}` gives stratum-specific conditional effects on the atom probability or odds scale, whereas a regression-style model such as

```math
\text{logit}\, p_r(l) = \alpha_0 + \alpha r + \alpha_l
```

imposes a common conditional log-odds treatment effect across strata.
- **Commentary / Preliminary Conclusion**: The key point is that stratifying by treatment arm and allowing categorical `L` to enter the mean is a valid and useful Bayesian construction, but it does not by itself yield a single adjusted conditional treatment effect. It yields one effect per stratum unless a no-interaction or common-shift restriction is added. Therefore there are really two different goals in this setting. If the scientific interest is in heterogeneous treatment effects across levels of `L`, then the more flexible arm-stratified formulation is preferable and the output should be stratum-specific posterior contrasts. If the goal is to reproduce the spirit of the frequentist adjusted model, then the Bayesian model should impose a common treatment effect across levels of `L` in both the atom component and the continuous component, with the DP mixture used only to model residual distributional flexibility around that common-shift structure. So the recommended interpretation is that “adjusted conditional effect” in the Bayesian DP framework requires a structural restriction, not just conditioning on the categorical covariate.

### Hierarchically Defined Model Over Levels of `L`

- **Name**: Hierarchical Bayesian model across categorical levels of `L`
- **Idea Description**: A natural compromise between full stratification and strict common-effect regression is to define a hierarchical model over the levels of a categorical covariate `L`. Suppose `L in {1, ..., K}`. For each treatment arm `r` and level `l`, write

```math
Y \mid R=r, L=l \sim p_{rl}\delta_a + (1-p_{rl})F^c_{rl}.
```

The key feature is that the cell-specific quantities `p_{rl}` and `F^c_{rl}` are not independent across `l`, but are tied together through higher-level priors. For the atom part one possibility is

```math
\text{logit}\, p_{rl} = \alpha_r + u_l^{(A)} + v_{rl}^{(A)},
```

with

```math
u_l^{(A)} \sim N(0,\sigma_A^2), \qquad v_{rl}^{(A)} \sim N(0,\tau_A^2),
```

so that levels of `L` share information through the common random effects. If one wants a common adjusted conditional treatment effect on the log-odds scale, one can instead write

```math
\text{logit}\, p_{rl} = \alpha_0 + \alpha r + u_l^{(A)},
```

and optionally add `v_{rl}^{(A)}` only if treatment heterogeneity across levels is to be allowed. For the continuous non-atom part, a parallel construction is

```math
Y_i \mid (Y_i \neq a, R_i=r, L_i=l, G_{rl})
\sim \int K(\cdot \mid m_{rl}, \theta)\, G_{rl}(d\theta),
```

with hierarchical mean structure such as

```math
m_{rl} = \beta_r + u_l^{(Y)} + v_{rl}^{(Y)},
```

or, under a common-shift restriction,

```math
m_{rl} = \beta_0 + \mu r + u_l^{(Y)}.
```

The distributional part can also be pooled hierarchically. A simple version is to let all `G_{rl}` share a common base measure with common hyperparameters. A richer version is a hierarchical Dirichlet process:

```math
G_r^\star \sim DP(\gamma_r H_r),
\qquad
G_{rl} \sim DP(\alpha_r G_r^\star),
```

so that the level-specific mixtures within arm `r` share atoms and borrow strength across levels of `L`. One could also define a single global `G^\star` above both arms if cross-arm borrowing is desired, though that is a stronger modeling choice.
- **Commentary / Preliminary Conclusion**: A hierarchically defined model over `L` seems especially attractive when `L` is categorical and some levels are sparse. It preserves the intuitive cell-specific interpretation of stratified models while avoiding the instability of fitting each `(R,L)` cell in isolation. It also lets us choose how much structure to impose. With only level-specific random intercepts `u_l`, the model yields a common adjusted conditional treatment effect and stays close to the frequentist adjusted estimand. With additional arm-by-level random terms `v_{rl}`, it becomes a model for treatment-effect heterogeneity across levels of `L`. For the continuous DP component, hierarchical pooling through shared hyperparameters or an HDP is a natural way to stabilize the level-specific mixtures without forcing them to be identical. This makes the hierarchical approach a strong candidate for the first categorical-covariate extension beyond simple stratification: it is more stable than separate cell-specific models, more interpretable than a fully predictor-dependent DP, and flexible enough to accommodate either common conditional effects or level-specific heterogeneity depending on the scientific goal. If `L` is ordinal rather than merely nominal, an additional refinement would be to replace exchangeable level effects by a smoothing prior, such as a random-walk prior over adjacent levels.

### Formal HDP Mixture Model for Categorical `L`

- **Name**: Formal hierarchical Dirichlet process mixture model for this problem
- **Idea Description**: Let `R \in \{0,1\}` denote treatment arm, let `L \in \{1, \ldots, K\}` be a categorical covariate, let `a` denote the atom, and let `C = S \setminus \{a\}` denote the non-atom support of the outcome space `(S, \mathcal{B})`. For each subject `i`, define `A_i = 1(Y_i \neq a)`. The full conditional response distribution in cell `(r,l)` is written as

```math
P(Y_i \in B \mid R_i=r, L_i=l)
= p_{rl}\,\delta_a(B) + (1-p_{rl})\,F_{rl}^c(B),
\qquad B \in \mathcal{B},
```

where `p_{rl} = P(Y=a \mid R=r, L=l)` and `F_{rl}^c` is a probability measure on `C`. The atom part can be modeled separately, for example through a hierarchical logistic prior such as

```math
\text{logit}\, p_{rl} = \alpha_0 + \alpha r + u_l^{(A)} + v_{rl}^{(A)},
```

with suitable exchangeable priors on `u_l^{(A)}` and `v_{rl}^{(A)}`. The hierarchical Dirichlet process is then used for the continuous non-atom component. Let `(Theta_c, \mathcal{T}_c)` be the parameter space of the non-atom kernel `K_c(\cdot \mid \theta)` on `C`, and let `H_0` be a base probability measure on `Theta_c`. For each treatment arm `r`, introduce an arm-level random probability measure

```math
G_r^\star \sim DP(\gamma_r H_0),
```

and for each level `l` of `L` within arm `r`, introduce a level-specific random probability measure

```math
G_{rl} \sim DP(\alpha_r G_r^\star).
```

The continuous non-atom distribution in cell `(r,l)` is then

```math
F_{rl}^c(B) = \int_{\Theta_c} K_c(B \mid \theta)\, G_{rl}(d\theta),
\qquad B \subseteq C,
```

so that, conditionally on `A_i=1`, `R_i=r`, and `L_i=l`,

```math
Y_i \mid G_{rl} \sim \int_{\Theta_c} K_c(\cdot \mid \theta)\, G_{rl}(d\theta).
```

If one wants a structured mean inside the non-atom component, the kernel can be centered at a cell-specific or constrained location parameter. For example,

```math
Y_i \mid (A_i=1, R_i=r, L_i=l, G_{rl})
\sim \int K_c(\cdot \mid m_{rl}, \theta)\, G_{rl}(d\theta),
```

where `m_{rl}` may be fully cell-specific or may satisfy a common-shift restriction such as `m_{rl} = \beta_0 + \mu r + u_l^{(Y)}`. A useful latent representation is obtained by introducing mixture atoms `\theta_{rlj}` and cluster labels `z_i`:

```math
G_r^\star = \sum_{h=1}^\infty \pi_{rh} \delta_{\phi_{rh}},
\qquad
G_{rl} = \sum_{h=1}^\infty \omega_{rlh} \delta_{\phi_{rh}},
```

where the cell-specific weights `\omega_{rlh}` are random and the atoms `\phi_{rh}` are shared within arm `r`. Then

```math
z_i \mid R_i=r, L_i=l \sim \sum_{h=1}^\infty \omega_{rlh}\,\delta_h,
\qquad
Y_i \mid (A_i=1, z_i=h, R_i=r, L_i=l) \sim K_c(\cdot \mid \phi_{rh}).
```

This representation makes clear that different levels of `L` within the same arm share mixture atoms through the common parent measure `G_r^\star`, which is the mechanism by which the HDP borrows strength across levels of `L`.
- **Commentary / Preliminary Conclusion**: In this problem, the HDP is most naturally applied only to the continuous non-atom component, while the atom probability remains in the explicit Bernoulli part of the model. This preserves the scientific interpretation of the atom as a distinct event and uses the HDP specifically to stabilize the non-atom distributions across categorical levels of `L`. Relative to independent cell-specific DPs, the main advantage is that sparse levels can reuse mixture components that are well supported in richer levels. The main modeling decision is whether the location structure `m_{rl}` is unconstrained or restricted: an unconstrained `m_{rl}` yields stratum-specific treatment effects, whereas a constrained common-shift specification can recover a single adjusted conditional effect analogous to the frequentist model. So the formal recommendation is that, when `L` is categorical and some levels are sparse, an explicit two-part model with an HDP prior on `F_{rl}^c` is a principled way to pool information across levels of `L` while preserving treatment-arm separation and the special role of the atom.

### Hierarchical Pooling Over `L` Versus Putting `L` in the Mean Structure

- **Name**: Compare hierarchical modeling across levels of `L` with regression on `L` in the mixture mean
- **Idea Description**: There are two conceptually different ways to use a categorical covariate `L` in the continuous non-atom part of the Bayesian model. One approach is hierarchical pooling over levels of `L`, where each cell `(R=r, L=l)` has its own conditional distribution or mixture `F_{rl}^c`, but those cell-specific distributions are linked through higher-level priors such as random effects or a hierarchical Dirichlet process. In that approach, `L` indexes related subpopulations, and the model borrows strength across those subpopulations while still allowing each of them to have its own distributional features. The other approach is to put `L` directly into the mean structure of the mixture model, for example through

```math
Y_i \mid (A_i=1, R_i=r, L_i=l, G_r)
\sim \int K_c(\cdot \mid m_r(l), \theta)\, G_r(d\theta),
```

with

```math
m_r(l) = \beta_{0r} + \beta_{lr}
```

or, under a common-shift restriction,

```math
m_r(l) = \beta_0 + \mu r + \beta_l.
```

In this second formulation, `L` is treated primarily as an explanatory variable that shifts the location of the non-atom distribution, while the residual mixture `G_r` captures deviations from a simple parametric error model. The first approach treats `L` as defining multiple related conditional distributions; the second treats `L` as entering a regression model for a single conditional location structure.
- **Commentary / Preliminary Conclusion**: From an interpretation perspective, hierarchical pooling over `L` is more natural when the scientific question is whether different levels of `L` correspond to genuinely different outcome distributions, possibly beyond simple shifts in mean. It preserves a cell-specific view of the data and makes it easy to talk about stratum-specific treatment effects and stratum-specific distributional differences. By contrast, including `L` in the mean structure is more natural when the goal is to adjust for `L` while preserving a lower-dimensional treatment effect, especially a common adjusted conditional effect analogous to the frequentist model. In that case, the role of `L` is mainly to explain baseline variation, not to define fundamentally distinct distributions. From an implementation perspective, regression on `L` in the mean structure is usually much simpler. It can often be fit with a single mixture model plus fixed or random effects, and it avoids introducing a separate cell-specific DP or HDP object for each level of `L`. Hierarchical pooling across `L` is more modular and often more faithful to the stratified scientific story, but it is computationally heavier, especially if the continuous part uses an HDP and the atom part also has hierarchical random effects. The hierarchical approach also creates more design choices: whether to pool only the means, only the mixture distributions, or both; whether to share information within arms only or across arms as well; and whether to allow treatment-effect heterogeneity across levels. So the tradeoff is fairly sharp. If the goal is stable estimation of related but potentially distinct stratum-specific distributions, hierarchical modeling over `L` is preferable. If the goal is a parsimonious adjusted conditional effect and straightforward computation, putting `L` into the mean structure is preferable. A useful practical distinction is therefore: use regression on `L` in the mean structure when adjustment is the main goal, and use hierarchical pooling over `L` when modeling heterogeneity across levels is itself scientifically important.

### Possible Implementations of the Mixture Construction

- **Name**: Practical implementation options for the DP or HDP mixture construction
- **Idea Description**: The mixture construction can be implemented in several equivalent computational forms. For the continuous non-atom component, the most direct representation is a latent-allocation mixture model:

```math
Y_i \mid (A_i=1, z_i=h) \sim K_c(\cdot \mid \theta_h),
\qquad
P(z_i=h \mid w) = w_h,
```

with the weights `w_h` and atoms `\theta_h` generated by a Dirichlet-process-type prior. One option is an infinite stick-breaking representation,

```math
v_h \sim Beta(1,\alpha), \qquad
w_h = v_h \prod_{g < h}(1-v_g),
\qquad
\theta_h \sim H,
```

possibly truncated at some large `H_max` for computation. This leads to a finite-but-flexible approximation that can be fit with standard blocked Gibbs or Hamiltonian Monte Carlo style methods if the truncation is explicit. A second option is a collapsed or partially collapsed Gibbs sampler based on the Chinese restaurant process or Pólya urn representation, where the weights are integrated out and the sampler updates cluster assignments directly. For ordinary DP mixtures this leads to Neal-style samplers, and for hierarchical Dirichlet processes the analogue is the Chinese restaurant franchise representation. A third option is slice sampling or retrospective sampling, which preserves the exact infinite representation while introducing auxiliary variables so that only finitely many components are active in each iteration. In the HDP case, the implementation may either follow the franchise representation or use a truncated stick-breaking construction at both the arm level and the level-within-arm level.

For this problem, the kernel `K_c` should be chosen so that it respects the support of the non-atom outcome. If the support is all of `\mathbb{R}`, Gaussian kernels are the easiest starting point. If the support is positive, log-normal or gamma kernels are more natural, or one may transform the data to `\mathbb{R}` and then use Gaussian kernels on the transformed scale. If the support is a closed interval, one can use transformed Gaussian kernels, beta kernels, or a latent Gaussian construction. If a structured mean `m_{rl}` is included, the implementation can center the kernel at that mean and place the DP or HDP on residual parameters such as local shifts, scales, or both. The atom part is best kept outside the mixture and handled separately through a Bernoulli or logistic component, so the mixture sampler only sees non-atom outcomes.
- **Commentary / Preliminary Conclusion**: For a package implementation, the most practical first choice is probably a truncated stick-breaking approximation with latent cluster labels, because it is explicit, stable, and easier to maintain than a fully nonparametric collapsed sampler. It also fits well with the two-part architecture: first sample the atom model, then sample the continuous non-atom mixture model conditional on the non-atom data. If the initial aim is the no-covariate model, a truncated DP mixture within each arm is likely enough. If the aim is hierarchical pooling across categorical levels of `L`, then a truncated HDP with arm-level and level-specific weights is a natural extension. Collapsed samplers and exact slice-based methods are theoretically elegant and may mix better in some settings, but they are harder to implement, debug, and expose reliably inside an R package. So the most realistic implementation path seems to be: start with support-appropriate kernels, latent allocations, and a truncated stick-breaking DP mixture for the continuous non-atom part; then extend to an HDP only when hierarchical pooling over levels of `L` becomes necessary. One additional design choice is whether to make the atoms `\theta_h` arm-specific or shared across arms. Sharing across arms can improve efficiency but weakens the separation between arm-specific outcome distributions, so for a first implementation it may be cleaner to keep separate arm-specific mixtures and only introduce shared atoms later if there is a clear reason to pool across arms.

### Estimands and Reporting Conventions

- **Name**: Core estimands for the Bayesian two-part model and how they should be reported
- **Idea Description**: The two-part Bayesian model naturally induces several distinct estimands, and the memo should distinguish them clearly. At the most basic level, for each treatment arm `r` and covariate value `L=l`, define the atom probability

```math
p_r(l) = P(Y=a \mid R=r, L=l),
```

the conditional non-atom mean

```math
\mu_r^c(l) = E(Y \mid Y \neq a, R=r, L=l),
```

and the conditional extended-outcome mean

```math
\mu_r^{ext}(l) = a\,p_r(l) + \{1-p_r(l)\}\mu_r^c(l).
```

These induce three natural conditional contrasts:

```math
\Delta_A(l) = p_1(l) - p_0(l),
```

```math
\Delta_c(l) = \mu_1^c(l) - \mu_0^c(l),
```

```math
\Delta_{ext}(l) = \mu_1^{ext}(l) - \mu_0^{ext}(l).
```

When the model imposes a common-shift or no-interaction restriction, some of these simplify to lower-dimensional estimands. For example, under

```math
\text{logit}\, p_r(l) = \alpha_0 + \alpha r + \alpha_l
```

the common conditional effect on the atom part is `\alpha` or `exp(\alpha)`, and under

```math
m_r(l) = \beta_0 + \mu r + \beta_l
```

the common conditional effect on the continuous non-atom mean is `\mu`. If no such restriction is imposed, then the correct estimands are stratum-specific quantities such as `\Delta_A(l)`, `\Delta_c(l)`, and `\Delta_{ext}(l)`. To obtain marginal estimands, the conditional quantities should be standardized over a reference distribution of `L`, for example the empirical distribution:

```math
\Delta_{ext}^{marg}
= E_L\{\Delta_{ext}(L)\}.
```

The same standardization can be applied to the atom and non-atom components separately. In a Bayesian analysis, each estimand should be summarized by posterior means or medians, credible intervals, and posterior probabilities such as `P(\Delta_{ext}^{marg} > 0 \mid data)` or `P(\Delta_{ext}(l) > 0 \mid data)`.
- **Commentary / Preliminary Conclusion**: The most important recommendation is that the memo should separate conditional, stratum-specific, and marginal estimands instead of treating “the treatment effect” as a single universal object. In the no-covariate model, this distinction collapses and the combined extended-outcome contrast is straightforward. Once covariates are present, however, the inferential target depends on the structure imposed by the model. If the model includes common-shift restrictions, it is coherent to report lower-dimensional adjusted conditional effects such as `\alpha` and `\mu`. If the model is stratified or hierarchical with treatment-effect heterogeneity, then the primary estimands should be level-specific posterior contrasts, possibly followed by a standardized marginal summary. For user-facing reporting, a sensible default is to report all three layers when available: atom effect, non-atom effect, and combined extended-outcome effect. The combined extended-outcome effect should probably be the headline estimand when a single summary is needed, while posterior probabilities and credible intervals should be reported alongside it rather than replaced by a p-value analogue.

### Model-Selection Guidance

- **Name**: Practical guidance for choosing among the Bayesian model variants
- **Idea Description**: The Bayesian framework discussed in this memo includes several nested modeling choices: whether to use no covariates or adjust for `L`, whether to treat `L` categorically or continuously, whether to aim for common adjusted conditional effects or heterogeneous stratum-specific effects, and whether to use independent cell-specific DPs, hierarchical pooling, or regression-plus-DP residual models. A practical model-selection rule should therefore depend on three things: the type of covariates, the scientific estimand of interest, and the amount of data available in each subgroup. In the simplest no-covariate case, the preferred model is separate treatment-arm explicit two-part modeling with a DP mixture for the continuous non-atom distribution in each arm. If `L` is categorical with a small number of well-populated levels and the scientific interest is in stratum-specific effects, then stratified explicit two-part models are the most transparent option. If `L` is categorical but some levels are sparse, hierarchical pooling over levels of `L` becomes preferable, with or without an HDP on the continuous component depending on how much distributional sharing is desired. If the scientific goal is instead to mimic the frequentist adjusted analysis and obtain a common conditional treatment effect, then the preferred model is one that puts `L` into the mean or linear predictor structure and imposes no interaction by default. If covariates are continuous or mixed, the first practical choice should usually be regression in both parts of the model, with a DP mixture on residuals in the continuous non-atom component, because this gives substantial flexibility without requiring a fully predictor-dependent DP. Fully covariate-dependent DPs or highly flexible hierarchical mixtures should be reserved for settings where the scientific question genuinely concerns distributional heterogeneity that cannot be captured by mean adjustment plus flexible residuals.
- **Commentary / Preliminary Conclusion**: A useful decision rule is the following. Use the no-covariate arm-specific two-part DP model when no adjustment is needed. Use stratified two-part models when `L` is purely categorical and the number of levels is small. Move to hierarchical pooling or an HDP when categorical levels are sparse and stratum-specific distributions remain scientifically important. Use regression on `L` in the atom and continuous mean structures when the primary goal is adjusted conditional treatment effects rather than detailed stratum-specific heterogeneity. For continuous or mixed covariates, prefer additive regression plus DP-mixture residuals as the first adjusted model, and only move to dependent-DP constructions if there is a clear need for richer conditional distributional effects. Finally, when sample size is limited or the inferential target is only a mean contrast, it may be better to use a simpler robust Bayesian model than to deploy the full DP or HDP machinery. So the model-selection guidance should emphasize that complexity is justified only when it serves a clear inferential purpose: the simplest model that matches the estimand and the data structure should generally be preferred.

### Complete Summary and Prioritized Recommendations

- **Name**: Overall summary of the memo and prioritized recommendations
- **Idea Description**: The memo has developed a Bayesian replacement for the current frequentist `TruncComp` logic by preserving the core two-part decomposition of the outcome into an atom indicator and a continuous non-atom component. The main conclusion is that the most coherent Bayesian starting point is an explicit two-part model, not a single undifferentiated model for the whole response. In the no-covariate setting, the recommended continuous model is a Dirichlet process mixture for the non-atom outcome within each treatment arm, while the atom probability is modeled separately. The memo also compared this explicit construction with an augmented-space DP formulation and concluded that the two can be formally equivalent under matched priors, but that the explicit parameterization is preferable in interpretation and implementation. When covariates are introduced, the memo distinguishes clearly between adjustment-oriented models and heterogeneity-oriented models. For categorical covariates, direct stratification is natural when cells are well populated, hierarchical pooling or an HDP becomes attractive when levels are sparse, and a common adjusted conditional effect requires a no-interaction or common-shift restriction. For continuous or mixed covariates, additive regression plus DP-mixture residuals is the preferred first extension because it preserves interpretable covariate effects while allowing flexible residual distributions. Across all variants, the memo emphasizes that the target estimand must be specified explicitly: atom effect, non-atom effect, combined extended-outcome effect, and whether these are conditional, stratum-specific, or marginal after standardization. Posterior probabilities and credible intervals, rather than Bayes factors by default, are the preferred reporting tools.
- **Commentary / Preliminary Conclusion**: The overall direction of the memo is now fairly coherent, and the main remaining task is to turn the conceptual structure into an implementation plan that does not overreach on the first version. The three main recommendations, in prioritized order, are:

1. Start with the simplest principled Bayesian model: an explicit two-part model with separate treatment-arm atom probabilities and arm-specific DP mixtures for the continuous non-atom distribution, using posterior summaries of the combined extended-outcome contrast as the main reported effect.
2. For covariate adjustment, choose the model according to the estimand rather than the other way around: use regression structures with common-shift restrictions when the goal is a frequentist-style adjusted conditional effect, and use stratified or hierarchical/HDP models only when heterogeneity across levels of `L` is scientifically important.
3. Keep the first implementation computationally conservative: use support-appropriate kernels, latent allocations, and a truncated stick-breaking DP mixture for the continuous non-atom part; postpone fully dependent DPs, highly customized exact samplers, and cross-arm sharing of mixture atoms until there is a clear inferential need.

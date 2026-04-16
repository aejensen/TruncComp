#' Internal interface normalization helpers
#'
#' Developer-facing documentation for the helpers that turn the public
#' `truncComp()` inputs into the standardized `Y`/`A`/`R` analysis frame used by
#' both estimation engines.
#'
#' @name trunccomp2-input-helpers
#' @title Internal interface normalization helpers
#' @description
#' These helpers validate the formula and default interfaces, normalize optional
#' covariate-adjustment inputs, infer the atom value when possible, and dispatch
#' into the shared fitting core.
#' @aliases normalizeFormulaAdjust
#' @aliases prepareFormulaAdjustment
#' @aliases prepareDefaultAdjustment
#' @aliases resolveDefaultAtom
#' @aliases truncComp_core
#' @usage
#' normalizeFormulaAdjust(adjust, outcome_name, treatment_name)
#' prepareFormulaAdjustment(data, adjust)
#' prepareDefaultAdjustment(adjust, n)
#' resolveDefaultAtom(y, a, atom = NULL)
#' truncComp_core(
#'   y, a, r, method, conf.level = 0.95, init = NULL,
#'   adjust_data = NULL, adjust_formula = NULL, atom = NULL
#' )
#' @details
#' ### `normalizeFormulaAdjust(adjust, outcome_name, treatment_name)`
#'
#' Validates a one-sided adjustment formula supplied through the public formula
#' interface. `adjust` must be `NULL` or a one-sided additive formula;
#' `outcome_name` and `treatment_name` are the variable names from the main
#' model formula and are used to reject collisions. The helper returns the
#' validated formula unchanged or `NULL` when no adjustment terms remain. It
#' errors on non-formulas, formulas using `.`, interaction terms, references to
#' the outcome or treatment variables, or the reserved standardized names `Y`,
#' `A`, and `R`. Its role is to ensure the later model-frame construction is
#' unambiguous.
#'
#' ### `prepareFormulaAdjustment(data, adjust)`
#'
#' Builds the adjustment-only model frame for the formula interface. `data` is
#' the user-supplied data frame and `adjust` is the validated one-sided formula.
#' The helper returns a list with `data` containing the covariate frame and
#' `formula` containing the original adjustment formula, or `NULL` components
#' when no adjustment variables are present. It fails if `na.fail` detects
#' missing values or if duplicate covariate names would make the later
#' standardized analysis data ambiguous. Its role is to keep adjustment handling
#' identical across `LRT` and `SPLRT`.
#'
#' ### `prepareDefaultAdjustment(adjust, n)`
#'
#' Normalizes adjustment covariates passed through `truncComp.default()`.
#' `adjust` may be `NULL`, a data frame, or a matrix, and `n` is the expected
#' number of rows matching `y`. The helper returns a list with standardized
#' adjustment data and an additive formula created with `reformulate()`. It
#' errors when row counts do not match, names are missing or duplicated, the
#' input is neither matrix nor data frame, or reserved names would collide with
#' the internal `Y`/`A`/`R` columns. Its role is to align the default interface
#' with the formula interface before both enter the same core fitter.
#'
#' ### `resolveDefaultAtom(y, a, atom = NULL)`
#'
#' Determines which atom value should be stored on the fitted object for the
#' default interface. `y` is the combined outcome vector, `a` is the observation
#' indicator, and `atom` is an optional explicit atom supplied by the caller.
#' The helper returns a single numeric atom. It errors when an explicit atom is
#' not scalar and finite, or when implicit inference is impossible because
#' `y[a == 0]` does not collapse to a single unique finite value. Its role is to
#' make later `Delta` calculations reproducible even when the user omits `atom`.
#'
#' ### `truncComp_core(y, a, r, method, conf.level = 0.95, init = NULL,
#' adjust_data = NULL, adjust_formula = NULL, atom = NULL)`
#'
#' Shared engine used by both public interfaces after argument normalization.
#' `y`, `a`, and `r` are the standardized outcome, observation, and treatment
#' vectors; `method` chooses `LRT` or `SPLRT`; `conf.level`, `init`, and `atom`
#' are stored on the result; and `adjust_data` plus `adjust_formula` represent
#' optional baseline covariates. The helper returns a `"TruncComp2"` object from
#' [LRT()] or [SPLRT()], augmented with the standardized analysis data and atom.
#' If `isDataOkay()` fails it returns a structured error object instead of
#' attempting estimation. Its role is to keep the two public interfaces thin and
#' to guarantee that both methods see the same canonical analysis frame.
#'
#' @seealso [truncComp()], [truncComp.default()], [LRT()], [SPLRT()]
#' @keywords internal
NULL

#' Internal object-construction helpers
#'
#' Developer-facing documentation for the low-level helpers that validate the
#' standardized data, summarize adjustment metadata, and create successful or
#' failed `"TruncComp2"` result objects.
#'
#' @name trunccomp2-object-helpers
#' @title Internal object-construction helpers
#' @description
#' These helpers centralize object creation and metadata normalization so every
#' estimation path returns the same top-level structure.
#' @aliases isDataOkay
#' @aliases adjustmentSpecification
#' @aliases isValid
#' @aliases truncCompMethod
#' @aliases newTruncComp2
#' @aliases returnErrorData
#' @usage
#' isDataOkay(d)
#' adjustmentSpecification(adjust)
#' isValid(truncCompObj)
#' truncCompMethod(method)
#' newTruncComp2(
#'   muDelta = NULL, muDeltaCI = NULL, alphaDelta = NULL, alphaDeltaCI = NULL,
#'   Delta = NULL, DeltaCI = NULL, DeltaMarginalCI = NULL,
#'   DeltaProjectedCI = NULL, DeltaProfileCI = NULL, W = NULL, p = NULL,
#'   method, conf.level, success, error = "", init = NULL, data = NULL,
#'   adjust = NULL, atom = NULL
#' )
#' returnErrorData(
#'   error, method, conf.level, init = NULL, data = NULL,
#'   adjust = NULL, atom = NULL
#' )
#' @details
#' ### `isDataOkay(d)`
#'
#' Performs the final minimum-data check on the standardized analysis frame `d`.
#' The helper inspects only the observed outcomes (`A == 1`) and returns `TRUE`
#' when both treatment groups contain at least two observed values, otherwise
#' `FALSE`. It does not throw; the caller decides whether to convert the result
#' into a warning or structured error object. Its role is to block likelihood
#' routines that require at least minimal within-group variation.
#'
#' ### `adjustmentSpecification(adjust)`
#'
#' Converts a validated adjustment formula into a compact human-readable string.
#' `adjust` is either `NULL` or an additive one-sided formula. The helper
#' returns `NULL` when there is no effective adjustment and otherwise collapses
#' the term labels with `" + "`. It is intentionally permissive once the formula
#' has passed earlier validation. Its role is to provide stable metadata for the
#' fitted object and printed summaries.
#'
#' ### `isValid(truncCompObj)`
#'
#' Small predicate for checking whether a result object represents a successful
#' fit. `truncCompObj` is any object expected to carry a `success` field. The
#' helper returns `TRUE` only when `success` is explicitly `TRUE`, and `FALSE`
#' otherwise. Its role is to give downstream helpers a single convention for
#' success checks.
#'
#' ### `truncCompMethod(method)`
#'
#' Normalizes short method identifiers into the descriptive labels stored on the
#' public object. `method` may be the short internal codes or already-expanded
#' labels. The helper returns the display string that should appear in summaries.
#' Unknown labels are passed through unchanged. Its role is to decouple the user
#' interface from internal branching strings.
#'
#' ### `newTruncComp2(...)`
#'
#' Constructs the canonical `"TruncComp2"` object. The arguments mirror the
#' fields exposed by the package: component estimates and intervals, `Delta`
#' derivatives, test statistics, method metadata, success state, optional
#' standardized data, and the fitted atom. The helper returns a named list with
#' class `c("TruncComp2", "list")`. It performs only light normalization,
#' mainly expanding `method` through `truncCompMethod()`. Its role is to ensure
#' the parametric and semi-parametric branches emit identically shaped results.
#'
#' ### `returnErrorData(error, method, conf.level, init = NULL, data = NULL,
#' adjust = NULL, atom = NULL)`
#'
#' Convenience wrapper for creating failed `"TruncComp2"` objects. `error`
#' supplies the message to store, while the remaining arguments preserve method
#' and input metadata. The helper returns the same class as a successful fit but
#' with `success = FALSE` and the estimate fields left `NULL`. Its role is to
#' make print, summary, and confidence-interval dispatch behave consistently
#' even after estimation fails.
#'
#' @seealso [newTruncComp2()], [summary.TruncComp2()], [print.TruncComp2()]
#' @keywords internal
NULL

#' Internal example-loading and simulation helpers
#'
#' Developer-facing documentation for the utilities that locate bundled example
#' files and build reproducible simulated datasets for tests, examples, and
#' validation scripts.
#'
#' @name trunccomp2-simulation-helpers
#' @title Internal example-loading and simulation helpers
#' @description
#' These helpers keep example data access stable across installed and source
#' layouts and provide compact wrappers around the public simulation interface.
#' @aliases .truncComp2ExtdataPath
#' @aliases .validateSimulationInputs
#' @aliases .drawObservedOutcome
#' @aliases .simulateTruncatedGroup
#' @aliases simTruncData
#' @usage
#' .truncComp2ExtdataPath(filename)
#' .validateSimulationInputs(n, f0, f1, pi0, pi1, atom = 0)
#' .drawObservedOutcome(generator, n, label)
#' .simulateTruncatedGroup(n, r, generator, probability, label, atom = 0)
#' simTruncData(
#'   n, mu0, mu1, pi0, pi1, sigma = 1, dist = "norm", df = 4, atom = 0
#' )
#' @details
#' ### `.truncComp2ExtdataPath(filename)`
#'
#' Resolves the installed path to a packaged example `.rds` file. `filename` is
#' the expected file name within `extdata`. The helper returns the first
#' matching path under either the installed-package layout or the source-like
#' `inst/extdata` layout. It errors if no candidate exists. Its role is to make
#' the public example loaders robust in both installed and development contexts.
#'
#' ### `.validateSimulationInputs(n, f0, f1, pi0, pi1, atom = 0)`
#'
#' Checks the user-supplied inputs for [simulateTruncatedData()]. `n` must be a
#' positive integer, `f0` and `f1` must be functions, `pi0` and `pi1` must lie
#' in `[0, 1]`, and `atom` must be scalar and finite. The helper returns `NULL`
#' invisibly and throws descriptive errors on invalid inputs. Its role is to
#' keep the public simulator and internal wrappers from repeating the same input
#' validation logic.
#'
#' ### `.drawObservedOutcome(generator, n, label)`
#'
#' Calls a group-specific outcome generator and standardizes its output.
#' `generator` is the user-provided function, `n` is the requested sample size,
#' and `label` names the generator in error messages. The helper returns a
#' numeric vector of length `n`. If the generator returns a single scalar, that
#' first scalar is kept and the helper obtains the remaining draws by repeated
#' calls to `generator(1)`. It errors when the result is not finite numeric data
#' of the correct length. Its role is to make the simulation interface flexible
#' while keeping downstream code simple.
#'
#' ### `.simulateTruncatedGroup(n, r, generator, probability, label, atom = 0)`
#'
#' Simulates one treatment arm. `n` is the group size, `r` is the treatment-arm
#' label to store in the `R` column, `generator` produces observed outcomes,
#' `probability` is the Bernoulli observation probability, `label` is passed to
#' `.drawObservedOutcome()`, and `atom` fills the missing outcomes. The helper
#' errors if an observed draw equals `atom`, because the package treats `atom`
#' as the unique code for an unobserved outcome. It returns a standardized data
#' frame with `R`, `A`, and `Y`. It inherits input validation failures from the
#' upstream helpers. Its role is to keep the public simulator symmetric across
#' the two treatment groups.
#'
#' ### `simTruncData(n, mu0, mu1, pi0, pi1, sigma = 1, dist = "norm", df = 4,
#' atom = 0)`
#'
#' Legacy convenience wrapper around [simulateTruncatedData()]. `mu0` and `mu1`
#' define the group means, `sigma` controls the Normal scale, `dist` chooses
#' either `"norm"` or the shifted `"t-sq"` generator, `df` supplies the degrees
#' of freedom for the `t` case, and the remaining arguments are passed through
#' to the public simulator. The wrapper validates that `mu0` and `mu1` are
#' scalar finite numeric values, that the Normal case uses a positive finite
#' `sigma`, and that the `t` case uses a finite `df > 1` so the documented group
#' means are well-defined. The helper returns the same standardized data frame
#' as [simulateTruncatedData()]. Its role is to preserve a compact simulation
#' shortcut used in older scripts and notes.
#'
#' @seealso [simulateTruncatedData()], [loadTruncComp2Example()],
#'   [loadTruncComp2AdjustedExample()]
#' @keywords internal
NULL

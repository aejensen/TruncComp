#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
manuscript_dir=$(cd "${script_dir}/.." && pwd)
source_file="${manuscript_dir}/manuscript.tex"
bundle_dir="${manuscript_dir}/self_contained"
document_stem="two_component_treatment_comparisons"

stage_root=$(mktemp -d "${manuscript_dir}/.self-contained-build.XXXXXX")
stage_root=$(cd "${stage_root}" && pwd -P)
stage_bundle="${stage_root}/self_contained"
previous_bundle="${stage_root}/previous_self_contained"
driver_dir="${stage_root}/drivers"
split_work_dir="${stage_root}/split_work"
verify_dir="${stage_root}/verify"
bibliography_build_dir="${stage_root}/bibliography_build"
bibliography_driver="${bibliography_build_dir}/manuscript_for_bibliography.tex"
bibliography_file="${bibliography_build_dir}/manuscript_for_bibliography.bbl"
supplement_flatten_input="${stage_root}/manuscript_for_supplement_bundle.tex"

inline_stem="${document_stem}_inline_proofs"
appendix_stem="${document_stem}_appendix_proofs"
main_stem="${document_stem}_main"
supplementary_stem="${document_stem}_supplementary"

cleanup() {
  if [[ -e "${previous_bundle}" && ! -e "${bundle_dir}" ]]; then
    mv "${previous_bundle}" "${bundle_dir}"
  fi

  if [[ "${KEEP_MANUSCRIPT_BUNDLE_TEMP:-0}" == "1" ]]; then
    printf 'Preserved temporary bundle workspace: %s\n' "${stage_root}" >&2
  else
    rm -rf "${stage_root}"
  fi
}
trap cleanup EXIT

for command_name in latexmk latexpand perl kpsewhich; do
  if ! command -v "${command_name}" >/dev/null 2>&1; then
    printf 'Required command not found: %s\n' "${command_name}" >&2
    exit 1
  fi
done

if [[ ! -f "${source_file}" ]]; then
  printf 'Canonical manuscript source not found: %s\n' "${source_file}" >&2
  exit 1
fi

proof_style_file=$(kpsewhich proof-at-the-end.sty)
if [[ -z "${proof_style_file}" || ! -f "${proof_style_file}" ]]; then
  printf 'The TeX installation does not provide proof-at-the-end.sty.\n' >&2
  exit 1
fi

flatten_driver() {
  local driver_file=$1
  local output_file=$2

  (
    cd "${manuscript_dir}"
    latexpand \
      --keep-comments \
      --fatal \
      --expand-bbl "${bibliography_file}" \
      --output "${output_file}" \
      "${driver_file}"
  )
}

prepare_flattened_source() {
  local flattened_file=$1

  perl -0pi -e '
    BEGIN { $count = 0; }
    $count += s{\\graphicspath\{\{build/figures/\}\}}{\\graphicspath{{figures/}}}g;
    END { die "Expected exactly one canonical graphicspath declaration\n" unless $count == 1; }
  ' "${flattened_file}"

  perl -0pi -e 's/[ \t]+$//mg; s/\s*\z/\n/' "${flattened_file}"

  if grep -Eq '\\(input|include|bibliography)[[:space:]]*\{' "${flattened_file}"; then
    printf 'Flattened source still contains an external LaTeX input or bibliography: %s\n' \
      "${flattened_file}" >&2
    grep -nE '\\(input|include|bibliography)[[:space:]]*\{' "${flattened_file}" >&2
    exit 1
  fi
}

copy_referenced_figures() {
  while IFS= read -r figure_path; do
    [[ -n "${figure_path}" ]] || continue

    case "${figure_path}" in
      /*|*..*)
        printf 'Unsafe figure path in flattened source: %s\n' "${figure_path}" >&2
        exit 1
        ;;
    esac

    source_figure="${manuscript_dir}/build/figures/${figure_path}"
    if [[ ! -f "${source_figure}" ]]; then
      printf 'Referenced figure not found: %s\n' "${source_figure}" >&2
      exit 1
    fi

    mkdir -p "${stage_bundle}/figures/$(dirname "${figure_path}")"
    cp "${source_figure}" "${stage_bundle}/figures/${figure_path}"
  done < <(
    perl -ne 'while (/\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}/g) { print "$1\n" }' \
      "$@" | sort -u
  )
}

write_driver() {
  local driver_file=$1
  shift

  printf '%s\n' "$@" > "${driver_file}"
  printf '\\input{%s}\n' "${source_file}" >> "${driver_file}"
}

printf 'Building the bibliography used by all four flattened sources...\n'
mkdir -p "${bibliography_build_dir}"
printf '\\input{%s}\n' "${source_file}" > "${bibliography_driver}"
(
  cd "${bibliography_build_dir}"
  TEXINPUTS="${manuscript_dir}//:" \
  BIBINPUTS="${manuscript_dir}//:" \
    latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
      "${bibliography_driver}"
)

if [[ ! -f "${bibliography_file}" ]]; then
  printf 'Bibliography build did not generate the expected file: %s\n' \
    "${bibliography_file}" >&2
  exit 1
fi

mkdir -p "${stage_bundle}/figures" "${driver_dir}"
cp "${proof_style_file}" "${stage_bundle}/proof-at-the-end.sty"
perl -pi -e 's/[ \t]+$//' "${stage_bundle}/proof-at-the-end.sty"

write_driver "${driver_dir}/${inline_stem}.tex" \
  '\newif\ifProofsToAppendix' \
  '\ProofsToAppendixfalse' \
  '\newif\ifProofsExternal' \
  '\ProofsExternalfalse'

write_driver "${driver_dir}/${appendix_stem}.tex" \
  '\newif\ifProofsToAppendix' \
  '\ProofsToAppendixtrue' \
  '\newif\ifProofsExternal' \
  '\ProofsExternalfalse'

write_driver "${driver_dir}/${main_stem}.tex" \
  '\newif\ifProofsToAppendix' \
  '\ProofsToAppendixtrue' \
  '\newif\ifProofsExternal' \
  '\ProofsExternaltrue' \
  '\newif\ifManuscriptIncludeSupplement' \
  '\ManuscriptIncludeSupplementfalse' \
  '\newif\ifManuscriptExternalMain' \
  '\ManuscriptExternalMaintrue'

printf 'Flattening the inline-proof, appendix-proof, and separate main sources...\n'
flatten_driver \
  "${driver_dir}/${inline_stem}.tex" \
  "${stage_bundle}/${inline_stem}.tex"
flatten_driver \
  "${driver_dir}/${appendix_stem}.tex" \
  "${stage_bundle}/${appendix_stem}.tex"
flatten_driver \
  "${driver_dir}/${main_stem}.tex" \
  "${stage_bundle}/${main_stem}.tex"

for flattened_file in \
  "${stage_bundle}/${inline_stem}.tex" \
  "${stage_bundle}/${appendix_stem}.tex" \
  "${stage_bundle}/${main_stem}.tex"; do
  prepare_flattened_source "${flattened_file}"
done

copy_referenced_figures \
  "${stage_bundle}/${inline_stem}.tex" \
  "${stage_bundle}/${appendix_stem}.tex" \
  "${stage_bundle}/${main_stem}.tex"

printf 'Generating the external proof stream with proof-at-the-end...\n'
cp -R "${stage_bundle}" "${split_work_dir}"
(
  cd "${split_work_dir}"
  latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
    "${main_stem}.tex"
)

external_proof_file="${split_work_dir}/${main_stem}-pratendmainproofs.tex"
if [[ ! -f "${external_proof_file}" ]]; then
  printf 'proof-at-the-end did not generate the expected proof stream: %s\n' \
    "${external_proof_file}" >&2
  exit 1
fi

external_proof_flatten_file="${stage_root}/${main_stem}-proofs-for-flattening.tex"
perl -0pe '
  s/\\begingroup\\renewcommand/\\begingroup\\def\\label#1{}\\renewcommand/g;
  s/\\makeatletter/\\makeatletter\n/g;
  s/\\makeatother/\\makeatother\n/g;
' "${external_proof_file}" > "${external_proof_flatten_file}"

PROOF_STREAM_PATH="${external_proof_flatten_file}" perl -pe '
  BEGIN {
    $proof_stream = $ENV{"PROOF_STREAM_PATH"};
    $proof_input = "\\input{" . $proof_stream . "}";
  }
  s/\\includeExternalAppendix\[mainproofs\]\{two_component_treatment_comparisons_main\}/$proof_input/g;
' "${source_file}" > "${supplement_flatten_input}"

{
  printf '%s\n' \
    '\newif\ifProofsToAppendix' \
    '\ProofsToAppendixtrue' \
    '\newif\ifProofsExternal' \
    '\ProofsExternaltrue' \
    '\newif\ifManuscriptIncludeMainText' \
    '\ManuscriptIncludeMainTextfalse' \
    '\newif\ifManuscriptExternalSupplement' \
    '\ManuscriptExternalSupplementtrue'
  printf '\\input{%s}\n' "${supplement_flatten_input}"
} > "${driver_dir}/${supplementary_stem}.tex"

printf 'Flattening the separate supplementary source and generated proof stream...\n'
flatten_driver \
  "${driver_dir}/${supplementary_stem}.tex" \
  "${stage_bundle}/${supplementary_stem}.tex"
prepare_flattened_source "${stage_bundle}/${supplementary_stem}.tex"
copy_referenced_figures "${stage_bundle}/${supplementary_stem}.tex"

printf 'Verifying all four variants with independent LaTeX builds...\n'
cp -R "${stage_bundle}" "${verify_dir}"
(
  cd "${verify_dir}"
  latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
    "${inline_stem}.tex"
  latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
    "${appendix_stem}.tex"
  latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
    "${main_stem}.tex"
  latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
    "${supplementary_stem}.tex"
  latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
    "${main_stem}.tex"
  latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
    "${supplementary_stem}.tex"
)

if grep -E 'LaTeX Warning: (Reference .* undefined|There were undefined references)|Package natbib Warning: (Citation .* undefined|There were undefined citations)' \
  "${verify_dir}"/*.log >/dev/null; then
  printf 'Undefined references or citations remain in a verified variant.\n' >&2
  grep -nE 'LaTeX Warning: (Reference .* undefined|There were undefined references)|Package natbib Warning: (Citation .* undefined|There were undefined citations)' \
    "${verify_dir}"/*.log >&2
  exit 1
fi

for verified_stem in \
  "${inline_stem}" \
  "${appendix_stem}" \
  "${main_stem}" \
  "${supplementary_stem}"; do
  cp "${verify_dir}/${verified_stem}.pdf" "${stage_bundle}/${verified_stem}.pdf"
done

unexpected_file=$(find "${stage_bundle}" -type f \
  ! -name '*.tex' ! -name '*.pdf' ! -name '*.sty' ! -name '*.png' \
  ! -name '*.jpg' ! -name '*.jpeg' ! -name '*.eps' -print -quit)
if [[ -n "${unexpected_file}" ]]; then
  printf 'Unexpected file in final bundle: %s\n' "${unexpected_file}" >&2
  exit 1
fi

if [[ -e "${bundle_dir}" ]]; then
  mv "${bundle_dir}" "${previous_bundle}"
fi

if ! mv "${stage_bundle}" "${bundle_dir}"; then
  if [[ -e "${previous_bundle}" ]]; then
    mv "${previous_bundle}" "${bundle_dir}"
  fi
  exit 1
fi

printf '%s\n' \
  'Built self-contained manuscript variants:' \
  "  ${bundle_dir}/${inline_stem}.tex" \
  "  ${bundle_dir}/${appendix_stem}.tex" \
  "  ${bundle_dir}/${main_stem}.tex" \
  "  ${bundle_dir}/${supplementary_stem}.tex" \
  'and the corresponding four PDF files.'

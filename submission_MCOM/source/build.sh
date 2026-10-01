#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
command -v pdflatex >/dev/null 2>&1 || { printf '%s\n' 'pdflatex is required. Install a TeX distribution with amsart, algorithm, algpseudocode, microtype, and hyperref.' >&2; exit 1; }
mkdir -p build
for name in main cover_letter reproducibility_index; do
  for pass in 1 2; do
    pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build "$name.tex" > "build/${name}_pass${pass}.txt" 2>&1 || { tail -60 "build/${name}_pass${pass}.txt" >&2; exit 1; }
  done
  if grep -Eq 'LaTeX Warning: (Reference|Citation).*undefined|There were undefined references|Overfull \\hbox|Missing character:' "build/$name.log"; then
    printf 'Review the warnings in build/%s.log\n' "$name" >&2
    exit 2
  fi
done
printf '%s\n' 'Built build/main.pdf, build/cover_letter.pdf, and build/reproducibility_index.pdf.'


dir="$(cd "$(dirname "$0")" && pwd)"

Rscript --verbose $dir/experimental_design.r $@ 2>&1

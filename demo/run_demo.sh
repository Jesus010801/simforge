#!/usr/bin/env bash
# SimForge talk demo — environment -> declarative spec -> guarded construction
# -> GROMACS validation -> provenance.  ~1 minute, no MD.
#
#   ./demo/run_demo.sh            # run it
#   ./demo/run_demo.sh --capture  # also write demo/expected_output.txt
set -euo pipefail
cd "$(dirname "$0")/.."

OUT="${SIMFORGE_DEMO_OUT:-demo_system}"
run() { echo; echo "\$ $*"; "$@"; }

rm -rf "$OUT"

run simforge doctor
echo; echo "----- the entire system, declaratively -----"
cat demo/system.yaml

run simforge build demo/system.yaml --out "$OUT" --verbose
run simforge validate-system "$OUT"
run simforge inspect-run "$OUT"

echo
echo "The two components each declared 'opls_800' for a different element;"
echo "SimForge namespaced the collision. See it in the provenance:"
run python -c "import json,sys; print(json.dumps(json.load(open('$OUT/provenance.json'))['atomtype_renames'], indent=2))"

if [[ "${1:-}" == "--capture" ]]; then
    # capture via a temp file so the redirect target does not itself make the
    # tree "dirty" while the demo runs
    _tmp="$(mktemp)"
    "$0" > "$_tmp" 2>&1 || true
    mv "$_tmp" demo/expected_output.txt
    echo "wrote demo/expected_output.txt"
fi

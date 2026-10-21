#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
#
# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Run the SHiP chain: geometry -> aegir (simulation) -> Shannon
# (digitisation) -> Trout (reconstruction). Every stage reads only the files
# the previous one wrote, so each can run in a different environment.
#
# Usage: scripts/run_chain.sh [output-dir]      (default: chain-output)
#
# Settings, from the environment:
#   EVENTS          events to simulate (default 20)
#   SEED            seed for the gun, Geant4 and digitisation (default 42)
#   POT             protons on target Shannon assumes per spill; its time
#                   offset is drawn from [0, 1.2e9 ns * POT / 4e13) (default
#                   1e7, i.e. 0 to 300 ns)
#   P_MIN, P_MAX    gun momentum range in GeV/c (default 10, 50)
#   MAX_THETA       gun opening angle in rad (default 0.01)
#   VERTEX_Z        gun vertex in mm, upstream of the upstream tagger (default 30000)
#   DETECTOR_IDS    how aegir is told the detector ids: positional (aegir
#                   <= 0.5) or explicit (later aegir) (default positional)
#   AEGIR_RUN, SHANNON_RUN, TROUT_RUN
#                   command prefix to run that stage's phlex with, e.g.
#                   "pixi run --manifest-path ../aegir/pixi.toml" (default: none,
#                   i.e. this environment)
#   AEGIR_WORKFLOWS, TROUT_WORKFLOWS
#                   the packages' workflow directories (default: the installed
#                   ones under $CONDA_PREFIX/share)
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
out=$(mkdir -p "${1:-chain-output}" && cd "${1:-chain-output}" && pwd)

: "${CONDA_PREFIX:?run inside the pixi environment (pixi run chain)}"
EVENTS=${EVENTS:-20}
SEED=${SEED:-42}
POT=${POT:-1e7}
P_MIN=${P_MIN:-10}
P_MAX=${P_MAX:-50}
MAX_THETA=${MAX_THETA:-0.01}
VERTEX_Z=${VERTEX_Z:-30000}
DETECTOR_IDS=${DETECTOR_IDS:-positional}
AEGIR_WORKFLOWS=${AEGIR_WORKFLOWS:-$CONDA_PREFIX/share/aegir/workflows}
TROUT_WORKFLOWS=${TROUT_WORKFLOWS:-$CONDA_PREFIX/share/trout/workflows}
field_lib=$SHIPFIELD_ROOT/share/field

# The other environments find the field maps and the geometry here too.
export SHIPFIELD_ROOT SHIPGEOMETRY_ROOT

stage() {
    local name=$1 runner=$2
    echo "--- $name"
    # Another environment brings its own plugins; don't let it find ours.
    local env_cmd=()
    [[ -z $runner ]] || env_cmd=(env -u PHLEX_PLUGIN_PATH)
    # The runner is a command prefix, so it is split on purpose.
    # shellcheck disable=SC2086
    (cd "$out" && "${env_cmd[@]}" $runner phlex -c "$name.json" >"$name.log" 2>&1) || {
        echo "$name failed, last lines of $out/$name.log:" >&2
        tail -n 30 "$out/$name.log" >&2
        exit 1
    }
}

echo "--- geometry"
(cd "$out" && build_geometry -o ship_geometry.db >geometry.log 2>&1)

jsonnet -J "$AEGIR_WORKFLOWS" -J "$field_lib" \
    --ext-str events="$EVENTS" --ext-str seed="$SEED" \
    --ext-str p_min="$P_MIN" --ext-str p_max="$P_MAX" \
    --ext-str max_theta="$MAX_THETA" --ext-str vertex_z="$VERTEX_Z" \
    --ext-str detector_ids="$DETECTOR_IDS" \
    --ext-str geometry="$out/ship_geometry.db" \
    --ext-str rntuple_file=sim.root --ext-str histo_file=sim_validation.root \
    "$here/workflows/sim.jsonnet" >"$out/sim.json"
stage sim "${AEGIR_RUN:-}"

jsonnet --ext-str spills="$EVENTS" --ext-str seed="$SEED" --ext-str pot="$POT" \
    --ext-str input_file=sim.root \
    --ext-str rntuple_file=digi.root --ext-str histo_file=digi_validation.root \
    "$here/workflows/digi.jsonnet" >"$out/digi.json"
stage digi "${SHANNON_RUN:-}"

# Shannon's output has no spill boundaries and Trout hands the whole file to
# every spill, so all events are reconstructed as one spill.
jsonnet -J "$TROUT_WORKFLOWS" -J "$field_lib" \
    --ext-str spills=1 --ext-str input_file=digi.root \
    --ext-str geometry="$out/ship_geometry.db" \
    --ext-str rntuple_file=reco.root --ext-str histo_file=reco_validation.root \
    "$here/workflows/reco.jsonnet" >"$out/reco.json"
stage reco "${TROUT_RUN:-}"

echo "--- outputs in $out"

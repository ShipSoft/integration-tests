# shellcheck shell=bash
# SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
#
# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Pixi activation script, sourced by `pixi run` / `pixi shell`.

# aegir, Shannon and Trout install their phlex plugins into $CONDA_PREFIX/lib.
export PHLEX_PLUGIN_PATH="${CONDA_PREFIX}/lib${PHLEX_PLUGIN_PATH:+:$PHLEX_PLUGIN_PATH}"

# Bare geometry-DB and field-map names resolve under these prefixes.
export SHIPGEOMETRY_ROOT="${SHIPGEOMETRY_ROOT:-$CONDA_PREFIX}"
export SHIPFIELD_ROOT="${SHIPFIELD_ROOT:-$CONDA_PREFIX}"

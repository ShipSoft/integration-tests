// SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// Stage 2 of the chain: Shannon digitises the hits of the aegir output.
// Each aegir event becomes one spill, so `spills` must equal the number of
// simulated events (read_sim_file reads entry i for spill i). `pot` sets the
// window of Shannon's per-spill time offset: 1.2e9 ns * pot / 4e13.
local ext(name) = std.extVar(name);

{
  stage: 'digitisation',
  driver: {
    cpp: 'generate_layers',
    layers: {
      spill: { parent: 'job', total: std.parseInt(ext('spills')) },
    },
  },

  sources: {
    rntuple_source: {
      cpp: 'read_sim_file',
      input_file: ext('input_file'),
      ntuple_name: 'events',
      layer: 'spill',
    },
  },

  modules: {
    digitise_hits: {
      cpp: 'digitise_hits',
      layer: 'spill',
      seed: std.parseInt(ext('seed')),
      pot: std.parseJson(ext('pot')),
    },
    digitised_output: {
      cpp: 'digitised_output_module',
      mode: 'digi',
      rntuple_file: ext('rntuple_file'),
      histo_file: ext('histo_file'),
      creator: 'digitise_hits',
      layer: 'spill',
    },
  },
}

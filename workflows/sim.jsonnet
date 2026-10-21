// SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// Stage 1 of the chain: aegir simulates a muon gun through the full SHiP
// geometry with the muon-shield and spectrometer field maps.
//
// Render with jsonnet (scripts/run_chain.sh does this), with aegir's
// workflow directory and $SHIPFIELD_ROOT/share/field on the import path:
//
//   jsonnet -J <aegir workflows> -J $SHIPFIELD_ROOT/share/field \
//     --ext-str events=20 --ext-str seed=1 ... workflows/sim.jsonnet
local field_maps = import 'field_maps.libsonnet';
local lib = import 'lib.libsonnet';

local ext(name) = std.extVar(name);

// Detector id per SHiP::detector_id, one pattern each. The straw-tube and
// timing-detector patterns match only the active volumes; the others match
// everything below the subsystem's top volume, envelopes included.
local detectors = [
  { pattern: '/SHiP/upstream_tagger/', detector: 'UpstreamTagger' },
  { pattern: '/SHiP/decay_volume/sbt/sensors', detector: 'SurroundTagger' },
  { pattern: '/SHiP/trackers/straw_gas', detector: 'StrawTubes' },
  { pattern: '/SHiP/calorimeter/', detector: 'Calorimeter' },
  { pattern: '/SHiP/timing_detector/bar', detector: 'TimingDetector' },
];

// aegir up to 0.5 numbers the detectors by their position in
// `sensitive_volumes`, so the list above is in SHiP::detector_id order.
// Later aegir takes the detector from each entry and rejects bare patterns
// it does not know. Select with --ext-str detector_ids=positional|explicit.
local sensitive_volumes =
  if ext('detector_ids') == 'positional' then [d.pattern for d in detectors]
  else if ext('detector_ids') == 'explicit' then detectors
  else error 'detector_ids must be positional or explicit';

{
  stage: lib.stage,
  driver: lib.driver(std.parseInt(ext('events'))),
  sources: {
    field: lib.field_map(field_maps.magnets),
    geometry: lib.geomodel_geometry {
      db_file: ext('geometry'),
      sensitive_volumes: sensitive_volumes,
    },
    gun: lib.gun {
      seed: std.parseInt(ext('seed')),
      p_min: std.parseJson(ext('p_min')),
      p_max: std.parseJson(ext('p_max')),
      max_theta: std.parseJson(ext('max_theta')),
      vertex_z: std.parseJson(ext('vertex_z')),
    },
  },
  modules: {
    geant4: lib.geant4 { seed: std.parseInt(ext('seed')) },
    // Keep events without hits: Shannon reads entry i for cell i, so every
    // generated event has to have an entry.
    output: lib.full_output(ext('rntuple_file'), ext('histo_file')) {
      filter_empty: false,
    },
  },
}

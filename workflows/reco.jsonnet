// SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
//
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// Stage 3 of the chain: Trout reconstructs the digitised hits, with the
// same geometry and spectrometer field map as the simulation. Needs Trout's
// workflow directory and $SHIPFIELD_ROOT/share/field on the import path.
local field_maps = import 'field_maps.libsonnet';
local reconstruction = import 'lib/reconstruction.libsonnet';

local ext(name) = std.extVar(name);
local spectrometer = [m for m in field_maps.magnets if m.name == 'Spectrometer'][0];

reconstruction({
  input_file: ext('input_file'),
  total_spills: std.parseInt(ext('spills')),
  db_file: ext('geometry'),
  field: {
    file: spectrometer.file,
    map: spectrometer.map,
    translation: spectrometer.translation,
  },
  rntuple_file: ext('rntuple_file'),
  histo_file: ext('histo_file'),
})

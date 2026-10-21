#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Check the outputs of scripts/run_chain.sh, stage by stage.

aegir:   one entry per simulated event, hits only with SHiP detector ids,
         and hits in the detectors a forward muon crosses.
Shannon: one collection per detector, each holding the hits aegir gave
         that detector, plus the simulation truth passed through unchanged.
Trout:   every reconstruction output present, the per-detector collections
         carrying Shannon's hits, the simulation truth passed through
         unchanged, and spectrometer tracks matching the generated muons
         (see check_tracks).

Prints every failure and exits non-zero if there was any. With
`--tracking report` the track checks are printed but do not fail the run.
"""

import argparse
import os
import sys

import numpy as np
import ROOT
from scipy.optimize import linear_sum_assignment

# SHiP::detector_id, and the collection each one ends up in at every stage.
DETECTORS = {
    0: ("UpstreamTagger", "ubt_hits", "upstream_tagger"),
    1: ("SurroundTagger", "sbt_hits", "surround_tagger"),
    2: ("StrawTubes", "straw_tubes_hits", None),
    3: ("Calorimeter", "calorimeter_hits", "calorimeter"),
    4: ("TimingDetector", "timing_detector_hits", "timing_detector"),
}
# Detectors every forward muon from the gun crosses. The surround tagger
# only sees secondaries, so it is not required to have hits.
REQUIRED = (0, 2, 3, 4)

ROOT.gInterpreter.Declare("""
ROOT::RVecD component(ROOT::RVec<ROOT::RVecD> const& v, int i) {
  ROOT::RVecD r;
  r.reserve(v.size());
  for (auto const& p : v) r.push_back(p[i]);
  return r;
}

double momentum_of(int track, ROOT::RVecI const& ids, ROOT::RVec<ROOT::RVecD> const& p) {
  for (std::size_t i = 0; i < ids.size(); ++i)
    if (ids[i] == track) return std::sqrt(p[i][0] * p[i][0] + p[i][1] * p[i][1] + p[i][2] * p[i][2]);
  return 0;
}
""")


class Checks:
    def __init__(self):
        self.errors = []

    def expect(self, ok, message):
        if not ok:
            self.errors.append(message)
        return ok


def entries(path, ntuple):
    return ROOT.RNTupleReader.Open(ntuple, path).GetNEntries()


def sim_hit_counts(path):
    """Hits per detector id in aegir's `events` RNTuple."""
    df = ROOT.RDataFrame("events", path).Define("ids", "ROOT::RVecI(sim_hits.detector_id)")
    ids = np.concatenate([np.asarray(v, dtype=int) for v in df.AsNumpy(["ids"])["ids"]] or [[]])
    values, counts = np.unique(ids, return_counts=True)
    return dict(zip(values.tolist(), counts.tolist()))


def check_sim(c, path, events):
    print(f"aegir: {path}")
    n = entries(path, "events")
    c.expect(n == events, f"aegir: {n} entries in 'events', expected {events}")
    counts = sim_hit_counts(path)
    for det_id, count in sorted(counts.items()):
        name = DETECTORS.get(det_id, ("not a SHiP detector",))[0]
        print(f"  detector {det_id} ({name}): {count} hits")
    unknown = sorted(set(counts) - set(DETECTORS))
    c.expect(not unknown, f"aegir: hits with detector ids {unknown}, Shannon only digitises 0-4")
    for det_id in REQUIRED:
        c.expect(counts.get(det_id, 0) > 0, f"aegir: no hits in {DETECTORS[det_id][0]}")
    n_particles = ROOT.RDataFrame("events", path).Define(
        "n", "sim_particles.track_id.size()").Sum("n").GetValue()
    return counts, int(n_particles)


# The fields of SHiP::SimHit, SHiP::SimParticle and SHiP::RecHit: scalars,
# then 3-vectors.
FIELDS = {
    "sim_hits": (["detector_id", "geometry_node_id", "track_id", "pdg_code",
                  "energy_deposit", "time", "path_length"],
                 ["position", "momentum"]),
    "sim_particles": (["track_id", "parent_id", "pdg_code", "energy", "time",
                       "creator_process"],
                      ["vertex", "endpoint", "momentum"]),
    "rec_hit": (["detector_id", "track_id", "pdg_code", "energy_deposit", "time",
                 "path_length"],
                ["position", "momentum"]),
}


def records(path, ntuple, kind, prefix):
    """Every record in `ntuple` as a row of its FIELDS[kind], rows sorted.

    `prefix` is the field holding the record: the collection name in aegir's
    `events`, which holds one vector of records per event, and for example
    "hit" or "hit.rec_hit" in the flat RNTuples, which hold one record per
    entry. Sorting makes comparisons independent of the record order.
    """
    scalars, vectors = FIELDS[kind]
    n_columns = len(scalars) + 3 * len(vectors)
    # RDataFrame (ROOT 6.40.04) crashes on an RNTuple without entries.
    if not entries(path, ntuple):
        return np.empty((0, n_columns))
    nested = ntuple == "events"
    df = ROOT.RDataFrame(ntuple, path)
    columns = [f"{prefix}.{name}" for name in scalars]
    for name in vectors:
        for i in range(3):
            df = df.Define(f"{name}_{i}", f"component({prefix}.{name}, {i})" if nested
                           else f"{prefix}.{name}[{i}]")
            columns.append(f"{name}_{i}")
    d = df.AsNumpy(columns)
    if nested:
        rows = np.column_stack([np.concatenate([np.asarray(v, dtype=float) for v in d[col]])
                                for col in columns])
    else:
        rows = np.column_stack([np.asarray(d[col], dtype=float) for col in columns])
    return rows[np.lexsort(rows.T[::-1])]


def compare_records(c, stage, what, source, read, read_expected):
    """Check that `read()` returns the records of `read_expected()`.

    A different number of records is left to the entry count checks.
    """
    try:
        got, expected = read(), read_expected()
    except Exception as e:  # noqa: BLE001 - ROOT raises plain exceptions
        c.expect(False, f"{stage}: cannot read {what}: {e}")
        return
    if got.shape != expected.shape:
        return
    differ = int(np.any(got != expected, axis=1).sum())
    c.expect(not differ, f"{stage}: {differ} of {len(got)} {what} differ from {source}")


def check_truth(c, stage, path, truth):
    """Compare the simulation truth `stage` passed through with aegir's."""
    for collection, expected in truth.items():
        compare_records(c, stage, collection, "aegir's",
                        lambda collection=collection: records(path, collection, collection,
                                                              "hit"),
                        lambda expected=expected: expected)


def check_digi(c, path, sim_counts, sim_particles, truth):
    print(f"Shannon: {path}")
    digi_counts = {}
    for det_id, (name, ntuple, _) in DETECTORS.items():
        try:
            n = entries(path, ntuple)
        except Exception as e:  # noqa: BLE001 - ROOT raises plain exceptions
            c.expect(False, f"Shannon: cannot open '{ntuple}': {e}")
            continue
        digi_counts[det_id] = n
        expected = sim_counts.get(det_id, 0)
        print(f"  {ntuple}: {n} hits from {expected} simulated")
        if expected:
            c.expect(0 < n <= expected,
                     f"Shannon: {n} hits in '{ntuple}' from {expected} simulated {name} hits")
        else:
            c.expect(n == 0, f"Shannon: {n} hits in '{ntuple}' without simulated {name} hits")
        if n:
            ids = ROOT.RDataFrame(ntuple, path).AsNumpy(["hit.rec_hit.detector_id"])
            wrong = set(np.unique(ids["hit.rec_hit.detector_id"]).tolist()) - {det_id}
            c.expect(not wrong, f"Shannon: '{ntuple}' holds hits of detector ids {sorted(wrong)}")
    n_hits = entries(path, "sim_hits")
    c.expect(n_hits == sum(sim_counts.values()),
             f"Shannon: {n_hits} sim_hits passed through, aegir wrote {sum(sim_counts.values())}")
    n_particles = entries(path, "sim_particles")
    c.expect(n_particles == sim_particles,
             f"Shannon: {n_particles} sim_particles passed through, aegir wrote {sim_particles}")
    check_truth(c, "Shannon", path, truth)
    return digi_counts


RECO_HISTOGRAMS = [
    "h_upstream_tagger_multiplicity",
    "h_surround_tagger_multiplicity",
    "h_spectrometer_track_multiplicity",
    "h_calorimeter_multiplicity",
    "h_timing_detector_multiplicity",
    "h_sim_hit_multiplicity",
    "h_sim_particle_multiplicity",
]


def check_reco(c, path, histo_path, digi_path, digi_counts, truth):
    print(f"Trout: {path}")
    for det_id, (name, digi_ntuple, ntuple) in DETECTORS.items():
        if ntuple is None:
            continue
        try:
            n = entries(path, ntuple)
        except Exception as e:  # noqa: BLE001
            c.expect(False, f"Trout: cannot open '{ntuple}': {e}")
            continue
        print(f"  {ntuple}: {n} hits")
        c.expect(n == digi_counts.get(det_id),
                 f"Trout: {n} hits in '{ntuple}', Shannon wrote {digi_counts.get(det_id)} "
                 f"to '{digi_ntuple}'")
        if n and n == digi_counts.get(det_id):
            compare_records(c, "Trout", f"hits in '{ntuple}'", f"Shannon's '{digi_ntuple}'",
                            lambda ntuple=ntuple: records(path, ntuple, "rec_hit",
                                                          "hit.rec_hit"),
                            lambda digi_ntuple=digi_ntuple: records(digi_path, digi_ntuple,
                                                                    "rec_hit", "hit.rec_hit"))
    for ntuple in ("sim_hits", "sim_particles"):
        n, expected = entries(path, ntuple), entries(digi_path, ntuple)
        c.expect(n == expected, f"Trout: {n} entries in '{ntuple}', Shannon wrote {expected}")
    check_truth(c, "Trout", path, truth)

    file = ROOT.TFile.Open(histo_path)
    if c.expect(file and not file.IsZombie(), f"Trout: cannot open {histo_path}"):
        for name in RECO_HISTOGRAMS:
            c.expect(bool(file.Get(name)), f"Trout: histogram '{name}' missing")
        file.Close()


# Track-matching thresholds, see check_tracks.
MATCH_TOLERANCE_MM = 1.0  # Shannon smears straw hits by 0.1 mm in x and y
# Centres in z of the four straw stations, each 1000 mm deep (centre
# ± 500 mm), as in the geometry's TrackersFactory.
STRAW_STATION_Z_MM = (84070.0, 86070.0, 93070.0, 95070.0)
STATION_HALF_DEPTH_MM = 500.0
MIN_STATIONS = len(STRAW_STATION_Z_MM)  # a reconstructable muon has hits in all of them
MIN_EFFICIENCY = 0.9
MAX_MOMENTUM_RESIDUAL = 0.1  # |p_reco / p_true - 1| of the track assigned to a muon
MIN_GOOD_MOMENTUM_FRACTION = 0.9


def truth_muons(path):
    """The gun muon of each aegir event: its momentum and straw-tube hits.

    Truth comes from aegir's output because Shannon's has no event
    boundaries. The gun muon is the event's primary, track_id 1.
    """
    df = (ROOT.RDataFrame("events", path)
          .Define("p_true", "momentum_of(1, sim_particles.track_id, sim_particles.momentum)")
          .Define("straw", "sim_hits.detector_id == 2 && sim_hits.track_id == 1")
          .Define("hx", "component(sim_hits.position, 0)[straw]")
          .Define("hy", "component(sim_hits.position, 1)[straw]")
          .Define("hz", "component(sim_hits.position, 2)[straw]"))
    d = df.AsNumpy(["p_true", "hx", "hy", "hz"])
    muons = []
    for p, x, y, z in zip(d["p_true"], d["hx"], d["hy"], d["hz"]):
        z = np.asarray(z)
        # Stations with at least one hit inside their envelope.
        stations = sum(bool(np.any(np.abs(z - centre) <= STATION_HALF_DEPTH_MM))
                       for centre in STRAW_STATION_Z_MM)
        muons.append({"p": float(p), "xy": np.column_stack([np.asarray(x), np.asarray(y)]),
                      "stations": stations})
    return muons


def assign_tracks(candidates):
    """Assign tracks to muons one to one.

    `candidates` holds a (chi2, track, muon) for every match. Of the
    assignments that give the most muons a track, takes the one with the
    lowest total chi2. Returns {muon: track}.
    """
    if not candidates:
        return {}
    muons = sorted({muon for _, _, muon in candidates})
    tracks = sorted({track for _, track, _ in candidates})
    row = {muon: i for i, muon in enumerate(muons)}
    column = {track: j for j, track in enumerate(tracks)}
    # A pair that does not match costs more than all matches together, so
    # the cheapest assignment is always one with the most matches.
    no_match = 1 + sum(chi2 for chi2, _, _ in candidates)
    cost = np.full((len(muons), len(tracks)), no_match)
    for chi2, track, muon in candidates:
        cost[row[muon], column[track]] = chi2
    rows, columns = linear_sum_assignment(cost)
    return {muons[i]: tracks[j] for i, j in zip(rows, columns) if cost[i, j] < no_match}


def check_tracks(c, reco_path, sim_path):
    """Match spectrometer tracks to the gun muons.

    Only tracks whose fit succeeded count. A track matches a muon if each of
    its measurements lies within MATCH_TOLERANCE_MM of one of the muon's
    straw-tube hits. Matches are then assigned one to one (see
    assign_tracks), so a track matching two muons counts for only one of
    them and does not keep another track from the other. Requires that at
    least MIN_EFFICIENCY of the muons crossing all MIN_STATIONS stations
    have a track, and that for at least MIN_GOOD_MOMENTUM_FRACTION of those
    the assigned track measures the momentum within MAX_MOMENTUM_RESIDUAL.
    """
    muons = truth_muons(sim_path)
    columns = ["hit.fit_status", "hit.q_over_p", "hit.chi2", "hit.input_measurements_x",
               "hit.input_measurements_y"]
    # RDataFrame (ROOT 6.40.04) crashes on an RNTuple without entries.
    if entries(reco_path, "spectrometer_tracks"):
        tracks = ROOT.RDataFrame("spectrometer_tracks", reco_path).AsNumpy(columns)
    else:
        tracks = {name: np.array([]) for name in columns}
    n_tracks = len(tracks["hit.q_over_p"])
    n_failed = int(np.sum(tracks["hit.fit_status"] != 0))
    reconstructable = [m for m in muons if m["stations"] >= MIN_STATIONS]
    print(f"  spectrometer_tracks: {n_tracks} tracks, {len(reconstructable)} of "
          f"{len(muons)} muons cross all {MIN_STATIONS} straw stations")

    candidates = []  # (chi2, track index, muon index) of every match
    for t, (status, chi2, xs, ys) in enumerate(zip(
            tracks["hit.fit_status"], tracks["hit.chi2"], tracks["hit.input_measurements_x"],
            tracks["hit.input_measurements_y"])):
        if status != 0:  # failed fit, no momentum
            continue
        meas = np.column_stack([np.asarray(xs), np.asarray(ys)])
        if not len(meas):  # would match any muon below
            continue
        for i, m in enumerate(reconstructable):
            dist = np.linalg.norm(meas[:, None, :] - m["xy"][None, :, :], axis=2)
            if len(m["xy"]) and np.all(dist.min(axis=1) < MATCH_TOLERANCE_MM):
                candidates.append((chi2, t, i))
    matched_tracks = len({t for _, t, _ in candidates})
    # muon index -> p_reco of its assigned track, q/p in 1/GeV
    best = {i: 1.0 / abs(tracks["hit.q_over_p"][t])
            for i, t in assign_tracks(candidates).items()}

    n = len(reconstructable)
    efficiency = len(best) / n if n else 0.0
    residuals = [abs(best[i] / reconstructable[i]["p"] - 1) for i in best]
    good = sum(r <= MAX_MOMENTUM_RESIDUAL for r in residuals)
    print(f"  tracking efficiency {efficiency:.2f} ({len(best)} of {n} muons), "
          f"{matched_tracks} matched ({len(best)} assigned), "
          f"{n_tracks - n_failed - matched_tracks} unmatched and "
          f"{n_failed} failed fits")
    if residuals:
        print(f"  momentum: {good} of {len(residuals)} assigned tracks within "
              f"{MAX_MOMENTUM_RESIDUAL:.0%}, median |p_reco/p_true - 1| "
              f"{np.median(residuals):.3f}")
    c.expect(n > 0, "Trout: no muon crossed all straw stations, check the gun")
    c.expect(efficiency >= MIN_EFFICIENCY,
             f"Trout: tracking efficiency {efficiency:.2f} below {MIN_EFFICIENCY}")
    c.expect(bool(best) and good >= MIN_GOOD_MOMENTUM_FRACTION * len(best),
             f"Trout: {good} of {len(best)} matched muons have the momentum within "
             f"{MAX_MOMENTUM_RESIDUAL:.0%}, need {MIN_GOOD_MOMENTUM_FRACTION:.0%}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_dir", nargs="?", default="chain-output")
    parser.add_argument("--events", type=int, default=int(os.environ.get("EVENTS", "20")))
    parser.add_argument("--tracking", choices=["require", "report"], default="require",
                        help="whether failed track checks fail the run (default: require)")
    args = parser.parse_args()

    def out(name):
        return os.path.join(args.output_dir, name)

    c = Checks()
    sim_counts, sim_particles = check_sim(c, out("sim.root"), args.events)
    truth = {collection: records(out("sim.root"), "events", collection, collection)
             for collection in ("sim_hits", "sim_particles")}
    digi_counts = check_digi(c, out("digi.root"), sim_counts, sim_particles, truth)
    check_reco(c, out("reco.root"), out("reco_validation.root"), out("digi.root"), digi_counts,
               truth)
    tracking = Checks()
    check_tracks(tracking, out("reco.root"), out("sim.root"))
    if args.tracking == "require":
        c.errors += tracking.errors
    else:
        for error in tracking.errors:
            print(f"KNOWN FAILURE (--tracking report) {error}")

    for error in c.errors:
        print(f"FAIL {error}", file=sys.stderr)
    print("all checks passed" if not c.errors else f"{len(c.errors)} check(s) failed")
    sys.stdout.flush()
    sys.stderr.flush()
    # Skip interpreter shutdown, which can crash in PyROOT's RNTuple teardown.
    os._exit(1 if c.errors else 0)


if __name__ == "__main__":
    main()

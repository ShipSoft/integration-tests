#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Check the outputs of scripts/run_chain.sh, stage by stage.

aegir:   one entry per simulated event, hits only with SHiP detector ids,
         and hits in the detectors a forward muon crosses.
Shannon: one collection per detector, each holding the hits aegir gave
         that detector, shifted in time by the per-spill offset, plus the
         simulation truth passed through.
Trout:   every reconstruction output present, the per-detector collections
         carrying Shannon's hits, and spectrometer tracks matching the
         generated muons (see check_tracks).

Prints every failure and exits non-zero if there was any. With
`--tracking report` the track checks are printed but do not fail the run.
"""

import argparse
import os
import sys

import numpy as np
import ROOT

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


# Shannon draws each spill's time offset from [0, window), with
# window = SPILL_TIME_NS * pot / NOMINAL_POT_PER_SPILL (digitise_hits.cpp).
SPILL_TIME_NS = 1.2e9
NOMINAL_POT_PER_SPILL = 4e13
# The mean offset of a collection must be at least this large. Shannon's
# default pot gives a window of only 0.3 ns.
MIN_MEAN_OFFSET_NS = 1.0


def mean_time(path, ntuple, column):
    """Mean of `column` over all entries, or None for an empty RNTuple."""
    if not entries(path, ntuple):
        return None
    return ROOT.RDataFrame(ntuple, path).Mean(column).GetValue()


def sim_mean_times(path):
    """Mean sim-hit time per detector id in aegir's `events` RNTuple."""
    df = (ROOT.RDataFrame("events", path)
          .Define("ids", "ROOT::RVecI(sim_hits.detector_id)")
          .Define("times", "ROOT::RVecD(sim_hits.time)"))
    d = df.AsNumpy(["ids", "times"])
    ids = np.concatenate([np.asarray(v, dtype=int) for v in d["ids"]] or [[]])
    times = np.concatenate([np.asarray(v, dtype=float) for v in d["times"]] or [[]])
    return {det_id: float(times[ids == det_id].mean()) for det_id in np.unique(ids).tolist()}


def check_time_offsets(c, digi_path, sim_path, pot):
    """Check that Shannon shifted each collection by its spill offsets.

    Shannon's output has no spill boundaries, so this compares means: the
    mean digitised time minus the mean simulated time of a detector is a
    weighted mean of the spill offsets, and must lie within the window.
    """
    window = SPILL_TIME_NS * pot / NOMINAL_POT_PER_SPILL
    sim_means = sim_mean_times(sim_path)
    for det_id, (_, ntuple, _) in DETECTORS.items():
        digi_mean = mean_time(digi_path, ntuple, "hit.rec_hit.time")
        if digi_mean is None or det_id not in sim_means:
            continue
        offset = digi_mean - sim_means[det_id]
        print(f"  {ntuple}: mean time offset {offset:.2f} ns (window {window:.1f} ns)")
        c.expect(MIN_MEAN_OFFSET_NS <= offset < window,
                 f"Shannon: mean time offset {offset:.3f} ns in '{ntuple}', expected "
                 f"{MIN_MEAN_OFFSET_NS} ns to {window:.1f} ns for pot {pot:g}")


def check_digi(c, path, sim_counts, sim_particles):
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


def check_reco(c, path, histo_path, digi_path, digi_counts):
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
    for ntuple in ("sim_hits", "sim_particles"):
        n, expected = entries(path, ntuple), entries(digi_path, ntuple)
        c.expect(n == expected, f"Trout: {n} entries in '{ntuple}', Shannon wrote {expected}")

    file = ROOT.TFile.Open(histo_path)
    if c.expect(file and not file.IsZombie(), f"Trout: cannot open {histo_path}"):
        for name in RECO_HISTOGRAMS:
            c.expect(bool(file.Get(name)), f"Trout: histogram '{name}' missing")
        file.Close()


# Track-matching thresholds, see check_tracks.
MATCH_TOLERANCE_MM = 1.0  # Shannon smears straw hits by 0.1 mm in x and y
MIN_STATIONS = 4  # straw stations a muon must cross to be reconstructable
MIN_EFFICIENCY = 0.9
MAX_MOMENTUM_RESIDUAL = 0.1  # |p_reco / p_true - 1| of the best track per muon
MIN_GOOD_MOMENTUM_FRACTION = 0.9
STATION_PITCH_MM = 2000.0  # the straw stations are 2 m apart in z


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
        # Station centres sit near multiples of the pitch (84070 mm, 86070 mm, ...)
        # and each station is only some 200 mm deep.
        stations = len(np.unique(np.round(z / STATION_PITCH_MM))) if len(z) else 0
        muons.append({"p": float(p), "xy": np.column_stack([np.asarray(x), np.asarray(y)]),
                      "stations": stations})
    return muons


def check_tracks(c, reco_path, sim_path):
    """Match spectrometer tracks to the gun muons.

    Only tracks whose fit succeeded count. A track belongs to a muon if each of its measurements lies within
    MATCH_TOLERANCE_MM of one of the muon's straw-tube hits. Requires that
    at least MIN_EFFICIENCY of the muons crossing all MIN_STATIONS stations
    have a track, and that for at least MIN_GOOD_MOMENTUM_FRACTION of those
    the best (lowest chi2) track measures the momentum within
    MAX_MOMENTUM_RESIDUAL.
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

    best = {}  # muon index -> (chi2, p_reco)
    matched_tracks = 0
    for status, qop, chi2, xs, ys in zip(*tracks.values()):
        if status != 0:  # failed fit, no momentum
            continue
        meas = np.column_stack([np.asarray(xs), np.asarray(ys)])
        for i, m in enumerate(reconstructable):
            dist = np.linalg.norm(meas[:, None, :] - m["xy"][None, :, :], axis=2)
            if len(m["xy"]) and np.all(dist.min(axis=1) < MATCH_TOLERANCE_MM):
                matched_tracks += 1
                p_reco = 1.0 / abs(qop)  # q/p in 1/GeV
                if i not in best or chi2 < best[i][0]:
                    best[i] = (chi2, p_reco)
                break

    n = len(reconstructable)
    efficiency = len(best) / n if n else 0.0
    residuals = [abs(best[i][1] / reconstructable[i]["p"] - 1) for i in best]
    good = sum(r <= MAX_MOMENTUM_RESIDUAL for r in residuals)
    print(f"  tracking efficiency {efficiency:.2f} ({len(best)} of {n} muons), "
          f"{matched_tracks} matched, {n_tracks - n_failed - matched_tracks} unmatched and "
          f"{n_failed} failed fits")
    if residuals:
        print(f"  momentum: {good} of {len(residuals)} best tracks within "
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
    parser.add_argument("--pot", type=float, default=float(os.environ.get("POT", "1e7")),
                        help="the pot Shannon ran with (default: $POT or 1e7)")
    parser.add_argument("--tracking", choices=["require", "report"], default="require",
                        help="whether failed track checks fail the run (default: require)")
    args = parser.parse_args()

    def out(name):
        return os.path.join(args.output_dir, name)

    c = Checks()
    sim_counts, sim_particles = check_sim(c, out("sim.root"), args.events)
    digi_counts = check_digi(c, out("digi.root"), sim_counts, sim_particles)
    check_time_offsets(c, out("digi.root"), out("sim.root"), args.pot)
    check_reco(c, out("reco.root"), out("reco_validation.root"), out("digi.root"), digi_counts)
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

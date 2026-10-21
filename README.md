# SHiP integration tests

Runs the SHiP software chain end to end and checks what each stage hands to
the next:

1. `build_geometry` (shipgeometry) writes the full SHiP geometry.
2. aegir simulates a muon gun through it, with the muon-shield and
   spectrometer field maps from ship-field-maps.
3. Shannon digitises the hits.
4. Trout reconstructs them.

Each package's own tests use synthetic input. These tests feed every stage
the real output of the stage before it, so they catch mismatches in file
layout, product names, detector ids and data-model versions between the
repositories.

## Running

```sh
pixi run test        # chain + validate
pixi run chain       # only run the chain, outputs in chain-output/
pixi run validate    # only check chain-output/
```

`scripts/run_chain.sh` reads its settings from the environment (`EVENTS`,
`SEED`, gun parameters and more; see the comment at its top). Each stage's
rendered configuration and log are kept in `chain-output/` next to its
output.

## What is checked

`scripts/validate_chain.py` checks, per stage:

- aegir writes one entry per event, only SHiP detector ids (0 to 4), and
  hits in the upstream tagger, straw tubes, calorimeter and timing detector.
- Shannon writes each detector's hits to that detector's collection, no
  more hits than aegir simulated, and passes the simulation truth through
  unchanged.
- Trout writes every reconstruction output and carries Shannon's hits into
  its per-detector collections.
- Trout's spectrometer tracks match the gun muons: at least 90% of the
  muons crossing all four straw stations have a track, and for at least 90%
  of those the best track's momentum is within 10% of the truth. Thresholds
  are constants at the top of `check_tracks`.

## Known problems

The track checks currently fail and run with `--tracking report`, which
prints them as known failures without failing the job:

- Trout 0.3.1 seeds only hit pairs whose times differ by less than 5 ns,
  but real hits in stations 1 and 2 are at least 6.7 ns apart (time of
  flight over 2 m), so no track is found.
- With the times corrected by hand, Trout finds every muon, but it pairs
  every station-1 hit with every station-2 hit and keeps the duplicates, and
  its placeholder 100 mm hit uncertainty leaves the momentum unmeasured.

Shannon's output has no spill boundaries and Trout hands the whole input
file to every spill, so the chain reconstructs all events as a single
spill.

Change the `validate` task to `--tracking require` once Trout reconstructs
real straw-tube hits.

## CI

`chain.yml` runs `pixi run test` on the released packages locked in
`pixi.lock`, on pull requests, on pushes to main and weekly.

`chain-mains.yml` builds aegir, Shannon and Trout from their main branches,
each in its own pixi environment, and runs the chain across them nightly.
The stages only exchange files, so they need not share an environment. This
job is advisory and requires the track checks, so it stays red until Trout
handles real hits.

aegir up to 0.5 numbers detectors by their position in `sensitive_volumes`;
later versions take the detector from each entry. `DETECTOR_IDS` selects
the form (`positional` for the released package, `explicit` for main).

## Contributing

Install the git hooks with `pixi run install-hooks`; `pixi run lint` runs
them on all files. Commits follow
[Conventional Commits](https://www.conventionalcommits.org). See
[AI_POLICY.md](AI_POLICY.md) for the use of AI tools.

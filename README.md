# oxDNA ANM regression tests

Regression tests for the protein models of oxDNA (the anisotropic network models, "ANM") and their
protein / nucleic-acid hybrids. They check that the current oxDNA reproduces the original
implementation in the `anm-oxdna` fork, in which these models were written, and they are meant to be
re-run during the refactor of the ANM code, then extended to CUDA and to oxNA.

| model | original code (`anm-oxdna`) | current oxDNA |
|---|---|---|
| protein ANM | `AC` | `AC` |
| protein ANMT (ANM + backbone bending/torsion) | `ACT` | `DNACT` with a protein-only topology (no `ACT`) |
| protein ANM + DNA (oxDNA2) | `DNANM` | `DNANM` |
| protein ANMT + DNA | `DNACT` | `DNACT` |
| protein ANM + RNA (oxRNA2) | `RNANM` | `RNANM` |
| protein ANMT + RNA | `RNACT` | `RNACT` |

## Running

```bash
# 1. build the two codes (see scripts/; both build out of tree and leave the checkouts untouched)
scripts/build_old.sh ~/software/anm-oxdna/oxDNA build/old             # -> build/old/bin/oxDNA
scripts/build_new.sh f802bea4 build/new                               # -> build/new/build/bin/oxDNA

# 2. run
./run_tests.py --old build/old/bin/oxDNA --new build/new/build/bin/oxDNA            # full CPU suite, ~14 min on an 8-core laptop
./run_tests.py --old build/old/bin/oxDNA --new build/new/build/bin/oxDNA --quick    # quick subset, ~30 s
./run_tests.py --old ... --new ... --only 'kdpg/*' 'ubq_dna/DNANM' --stages energy0,det
./run_tests.py --old ... --new ... --list                                           # list the cases
```

`--new` can be any oxDNA binary (e.g. a build of the refactoring branch). Results are written to
`results/<backend>[_quick]/` (`--out` to change): `summary.md` (one row per check), `results.json`
(all numbers, including per-residue RMSFs and per-frame trajectory deviations) and one run directory
per simulation (`<system>/<case>/<stage>/<old|new>[/repN]/` with the generated `input`, `log.dat`,
`energy.dat`, `trajectory.dat`). The exit status is non-zero if any check FAILs or ERRORs.
Check statuses: `PASS`, `FAIL`, `ERROR` (a run crashed or the analysis failed), `XFAIL`/`XPASS`
(failure expected / not expected, see `xfail` in the configs), `SKIP`.

`--jobs` (default: all cores) sets how many simulations run at once. CUDA: see `cuda_tests.md`.

`reference_results/CPU_f802bea4.md` is the summary of the full CPU suite for oxDNA f802bea4 against
anm-oxdna 015a855 (all checks pass except `kdpg/AC`, see "Known differences").

## What is tested

Each system lives in `systems/<name>/` with its inputs and a `config.json` that lists the cases
(model x system), the inputs given to each code, the run lengths and the expected failures.

| system | contents | cases |
|---|---|---|
| `ubiquitin` | PDB 1UBQ, 76 residues, 1 chain, 13 A network | AC, ANMT (ACT vs DNACT), DNACT, DNANM, RNANM, RNACT (protein only) |
| `cln025` | PDB 5AWL, 10-residue mini-protein, sparse 8 A network, per-spring kb = kt = 30 in the par file | AC, DNACT, RNACT (global kb/kt: the original RNACT cannot read per-spring constants) |
| `kdpg` | KDPG trimer from the original ANMUtils example, 605 residues, 3 chains | AC, ANMT (ACT vs DNACT), DNANM (protein only) |
| `ubq_dna` | ubiquitin + 12-bp DNA duplex held in contact by a mutual trap | DNANM, DNACT |
| `ubq_rna` | ubiquitin + 12-bp A-form RNA duplex held in contact by a mutual trap | RNANM, RNACT, each also with Debye-Hueckel switched off |
| `cage` | the hybrid cage example (KDPG + 393-nt DNA cage, mutual traps) | DNANM, DNACT |

All runs: MD, `thermostat = john`, T = 300 K, dt = 0.002, double precision on the CPU.

For every case the runner performs four checks:

1. **inputs** - the original code takes the protein bonds from the neighbour lists in the topology
   (plus the C-side neighbour), the current one from the parameter file: every spring of the `.par`
   file must also be a topology bond, otherwise the two codes simulate different networks.
2. **energy0** - potential energy (per particle) of a randomly perturbed configuration
   (`perturbed.dat`; the generated configurations sit exactly in the network minimum, where every
   term is ~0), from a `steps = 0` run of each code. Tolerance 2e-6 per particle (the energy file has
   6 decimals). The protein terms (springs, protein excluded volume, bending/torsion) are also
   evaluated independently with numpy (`tools/anm_reference.py`); for protein-only systems the engine
   energies must match it within the same tolerance, for hybrids the remainder (nucleic-acid +
   protein-nucleic terms) is reported. For RNA hybrids the step-0 energy is also computed with the
   Debye-Hueckel term switched off (`dh_strength = 1e-12`), which isolates the DH energy of each code.
3. **det** - a short run with the same seed in both codes. The two codes consume random numbers in the
   same order, so trajectories are identical up to round-off; the maximum position deviation of every
   printed frame must stay below `tol`. This is the strongest test; see "Tolerances" for why it is
   looser for nucleic acids.
4. **stat** - independent replicas (different seeds; 4-8 per code) started from an equilibrated
   configuration (`equilibrated.dat`, written by `tools/equilibrate.py`). Compared: the mean potential
   energy and the per-residue RMSF of the protein (Kabsch-aligned on the protein, Angstrom; 1 length
   unit = 8.518 A). Errors are the standard error of the replica means. Fails if
   |dE| > z_energy * sigma (z_energy = 4, 4.5 for nucleic systems), or if more than
   max(1, 5 %) residues have |dRMSF| > 4 sigma, or <z^2> over residues > 4, or the mean RMSF differs by
   more than 4 sigma + 2 %. The nucleic-acid RMSF is reported but not tested (fraying of the duplex
   ends is rare and dominates it). For protein-only systems the simulated RMSF is also compared with
   the ANM prediction of the parameter generator (`anm_msf.dat`; informational: correlation and ratio).

`--quick` runs `ubiquitin/{AC,ANMT,DNANM}`, `cln025/DNACT`, `ubq_dna/DNANM` and `ubq_rna/RNANM`
with shorter `det`/`stat` runs.

### Tolerances and why

- energy0: 2e-6 per particle = round-off of the 6-decimal energy output. On CUDA 1e-4 relative.
- det, protein only: 1e-6 (simulation units). ANM cases agree to ~1e-14. ANMT cases agree to
  ~1e-8: the original `LR_vector::normalize()` uses `sqrtf` even in double precision
  (`anm-oxdna/oxDNA/src/defs.h`), so the unit vector used by the bending term carries ~1e-7 relative
  error (switching off bending, but not torsion, makes the ANMT trajectories agree to 1e-15).
- det, nucleic acids: 2e-4 over 500 steps (cage: 2e-3 over 300 steps). The same `sqrtf` seeds 1e-7
  differences in the oxDNA/oxRNA terms, and nucleic-acid dynamics amplify them chaotically
  (~1e-5 after 500 steps for the duplexes, ~5e-4 for the cage).
- stat: statistical, see above. Replicas instead of block averages because the slow processes of the
  hybrids (fraying, protein-duplex contacts) make block errors too optimistic.

## Known differences between the original and the current code

1. **AC with several chains (current-code bug).** `ANMInteraction::pair_interaction_bonded` computes
   `_computed_r = q->pos - p->pos` without the minimum image. oxDNA shifts every strand's centre of
   mass into the box when it reads a configuration (and periodically with `fix_diffusion`), so for a
   protein made of several chains the inter-chain springs see distances off by a box vector.
   `kdpg/AC`: step-0 energy 1720.74 per particle (original and numpy reference: 0.0821), the run
   blows up to NaN. Single-chain proteins and DNANM/DNACT (which use `min_image`) are not affected;
   `fix_diffusion = 0` hides it. The `test/DNANM/new_topology/kdpg_AC` reference value in oxDNA
   (2164.606502) is this wrong number. Left as a FAIL in this suite.
2. **Protein torque from protein-nucleotide excluded volume (DNACT/RNACT).** The original adds the
   nucleotide's excluded-volume torque to the protein particle as well
   (`protein->torque -= protein->orientationT * -torquenuc`); the current code adds
   `(-r_hat * sigma/2) x F` about the protein centre, which is identically zero because the force is
   central (and uses `p->orientationT` instead of `protein->orientationT`, harmless only because the
   term is zero). The current behaviour conserves angular momentum, the original does not. Effect:
   identical energies, ANMT trajectories diverge at the first protein-nucleotide contact
   (`ubq_dna/DNACT`, `cage/DNACT` det are XFAIL; the deviation starts on the residue in contact and
   its backbone neighbours); averages (stat) agree within error.
3. **Debye-Hueckel in RNANM/RNACT.** The double counting of the DH term for RNA-RNA pairs was in the
   *pre-refactor current* oxDNA (`RNANMInteraction::pair_interaction_nonbonded` added `_debye_huckel`
   on top of `RNA2Interaction::pair_interaction_nonbonded`, which already includes it; removed in
   oxDNA commit 1f9f4e75). The original `anm-oxdna` RNANM/RNACT count it once. The `ubq_rna` cases
   measure it: the DH energy (E - E with DH off) is identical in both codes (ratio 1.000), and the
   total energies, trajectories and averages agree.
4. **Precision of the original ANMT.** See the det tolerances above (`sqrtf` in `normalize`).
5. **Input formats.** The original needs the protein neighbour lists in the topology; the current
   code ignores them and uses the `.par` file (hence the `inputs` check). The original RNACT and ACT
   only accept global `bending_k`/`torsion_k`; DNACT (both codes) and the current RNACT also accept
   per-spring kb kt in the par file. The original DNACT rejects par-file indices larger than the
   number of amino acids, so in hybrid topologies the protein must come first.

## Layout

```
run_tests.py            the runner
regtest/                file readers and statistics (block averages, Kabsch alignment, RMSF)
systems/<name>/         inputs + config.json for each system
  conf.dat              starting (network-minimum) configuration
  perturbed.dat         energy0 configuration (tools/perturb_conf.py)
  equilibrated.dat      start of the stat replicas (tools/equilibrate.py)
  anm.par / anmt.par    ANM and ANMT parameter files
  anm.top / mixed.top / hybrid.top   old-style topologies with neighbour lists
  anm_msf.dat           ANM-predicted fluctuations (protein-only systems)
tools/                  generators (see below)
scripts/                build scripts
cuda_tests.md           instructions for the GPU runs
```

### How the inputs were made

```bash
python tools/make_protein_system.py systems/pdb_downloads/1UBQ.pdb systems/ubiquitin
python tools/make_protein_system.py systems/pdb_downloads/5AWL.pdb systems/cln025 --cutoff 8 --kb 30 --kt 30
awk 'NR==1{print;next} NF==11{NF=9} {print}' systems/cln025/anmt.par > systems/cln025/anmt_global.par
# kdpg: kdpgt.dat, kdpg.par, kdpgt.par, kdpg.top from anm-oxdna/oxDNA/ANMUtils/examples/KDPG
# cage: acc2.top, ac.par, acc2_relaxed.dat, acforceponly.txt from oxDNA/examples/DNANM/hybrid_cage, plus
python tools/add_anmt_angles.py systems/cage/conf.dat systems/cage/anm.par systems/cage/anmt.par
python oxDNA/utils/generate-sa.py 20 seq.txt          # seq.txt: "DOUBLE GATTACAGCTCG" -> systems/ubq_dna/dna_duplex.*
python tools/make_rna_duplex.py GGAUCCAGUCAG systems/ubq_rna/rna_duplex
python tools/make_hybrid.py systems/ubiquitin systems/ubq_dna/dna_duplex.dat systems/ubq_dna/dna_duplex.top systems/ubq_dna
python tools/make_hybrid.py systems/ubiquitin systems/ubq_rna/rna_duplex.dat systems/ubq_rna/rna_duplex.top systems/ubq_rna
python tools/perturb_conf.py systems/X/conf.dat systems/X/perturbed.dat      # cln025: --sigma 0.03 --sigma-rot 0.2
python tools/equilibrate.py <current oxDNA binary> ubiquitin cln025 kdpg ubq_dna ubq_rna cage
```

`tools/models.py` is a copy of `ANMUtils/models.py` with three changes, all marked `PATCHED`:
`float(np.ravel(1 / slope)[0])` (numpy >= 2 refuses `float()` of a 1-element array), a headless
matplotlib backend, and normalised a3 vectors (the original writes `a3 = a1 x e_x` without
normalising it, while oxDNA normalises orientations on input, so the torsion reference
`d0 = a3_i . a3_j` in the par file did not correspond to the simulated configuration).

## Adding tests

- A new case: add an entry to `cases` in a system's `config.json`. Keys: `old`, `new` and `common`
  (input-file keys for each code / both; relative file names are resolved in the system directory),
  `reference` (par file and kb/kt for the numpy check), `quick`, `xfail`
  (`{"<check>": reason, "<check>@CUDA": reason, "*@CUDA": reason}`), `energy0_variants` /
  `stat_variants` (extra runs with modified inputs).
- A new system: a directory with the input files and a `config.json` (copy one). `protein_range` and
  `nucleic_range` select the particles for the RMSF; `cost` (ms per step) only orders the jobs.
- **oxNA / CPU vs CUDA of the same code.** Run the same binary family on both sides with
  `--old-inputs new` so both sides use the `new` inputs, e.g.
  `./run_tests.py --old oxDNA_cpu --new oxDNA_cuda --backend CUDA --old-inputs new --only 'oxna_*'`.
  The `det` stage is skipped across backends; `energy0` and `stat` (mean energy and RMSF, which also
  works for nucleic-acid-only systems by setting `protein_range` to the particles to analyse) are
  what such a comparison needs. A system for the `NA` interaction only needs its inputs and a
  `config.json` whose cases set `interaction_type = NA` in `common`.

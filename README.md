# oxDNA ANM regression tests

Regression tests for the protein models of oxDNA (the anisotropic network models, "ANM"), their
protein / nucleic-acid hybrids, and oxNA (`NA`). On the CPU they check that the current oxDNA
reproduces the original implementation in the `anm-oxdna` fork, where these models were written. For
the oxNA-based models (`NANM`/`NACT`, `NA`), which have no original implementation, they check
equivalence with the corresponding current models or compare with stored reference values. On the GPU
they compare the CUDA build with the CPU build of the same code (see `cuda_tests.md`).

| model | original code (`anm-oxdna`) | current oxDNA |
|---|---|---|
| protein ANM | `AC` | `AC` |
| protein ANMT (ANM + backbone bending/torsion) | `ACT` | `ACT` (also `DNACT` with a protein-only topology) |
| protein ANM + DNA (oxDNA2) | `DNANM` | `DNANM` |
| protein ANMT + DNA | `DNACT` | `DNACT` |
| protein ANM + RNA (oxRNA2) | `RNANM` | `RNANM` |
| protein ANMT + RNA | `RNACT` | `RNACT` |
| protein ANM/ANMT + DNA, RNA and DNA/RNA hybrids (oxNA) | - | `NANM`, `NACT` |
| DNA + RNA + hybrids (oxNA) | - | `NA` |

## Running

```bash
# 1. build the two codes (see scripts/; both build out of tree and leave the checkouts untouched)
scripts/build_old.sh ~/software/anm-oxdna/oxDNA build/old             # -> build/old/bin/oxDNA
scripts/build_new.sh anm_refactor build/new                           # -> build/new/build/bin/oxDNA

# 2. run
./run_tests.py --old build/old/bin/oxDNA --new build/new/build/bin/oxDNA            # full CPU suite, ~12 min on an 8-core laptop
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

`--jobs` (default: all cores) sets how many simulations run at once. CUDA (`--same-code
--backend CUDA --cuda-precision mixed|double`, `--smoke`, `scripts/cluster_cuda.sh`): see `cuda_tests.md`.

Reference results (full CPU suite against anm-oxdna 015a855):
- `reference_results/CPU_ca0afa94.md`: oxDNA anm_refactor ca0afa94. Everything passes; the only XFAILs
  are the two DNACT hybrid `det` checks (difference 2 below).
- `reference_results/CPU_f802bea4.md`: f802bea4. `kdpg/AC` fails (difference 1).

### Kinds of cases

- `"old_code": "original"` (default): `old` inputs with the `--old` binary (the original code) vs `new`
  inputs with `--new`.
- `"old_code": "current"`: both sides run with the `--new` binary, e.g. `NANM` with only DNA strands vs
  `DNANM` (`ubq_dna_nt`, `ubq_rna_nt`). These must agree bit for bit (`energy0_tol = 0`, `det` tol = 0).
- `"old_code": "stored"`: no reference code (`ubq_hybrid` NANM/NACT, the `na_*` oxNA systems). The test
  side is compared with `systems/<s>/stored.json`, written by `--update-stored` (energy0 value, replica
  averages and RMSFs, NVE drift). The stored values were made with seeds 1000+ and the checks run with
  2000+, so `stat` is a genuine statistical comparison. Regenerate them only on purpose, e.g. after an
  intended model change: `./run_tests.py --old X --new NEW --only 'ubq_hybrid/*' 'na_*' --update-stored`.
- `--same-code`: both sides use the `new` inputs, plus `new@CUDA` substitutions such as AC -> DNANM. The
  reference is a CPU build of the same code. Used for CPU vs CUDA.

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
| `ubq_dna_nt`, `ubq_rna_nt` | `ubq_dna`/`ubq_rna` with a new-style topology (`tools/to_new_topology.py`) | NANM vs DNANM/RNANM, NACT vs DNACT/RNACT (bit for bit) |
| `ubq_hybrid` | ubiquitin + DNA:RNA hybrid duplex with a DNA overhang (oxDNA `test/CUDA/NA/hybrid_overhang`), mutual trap | NANM, NACT (stored) |
| `na_dna_duplex`, `na_rna_duplex`, `na_hybrid_overhang`, `na_mixed` | the oxNA systems of oxDNA `test/CUDA/NA` (22-179 nt), T = 25 C, 0.5 M | NA (stored; CPU vs CUDA) |

`ubiquitin` and `kdpg` also have an `ACT` case: original ACT vs current ACT.

All runs: MD, `thermostat = john`, T = 300 K, dt = 0.002, double precision on the CPU.

For every case the runner performs up to five checks (stages):

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
5. **nve** - 20 000 steps at dt = 0.001 with `thermostat = no`, from `equilibrated.dat`. The conserved
   energy is U + K + ext/2: `force_energy` counts each mutual-trap pair twice. Its drift (linear fit
   over the run) and its fluctuations about the fit must satisfy test <= 3 x ref + 2e-4 per particle.
   On CUDA this directly tests the GPU forces, because the energies are computed on the CPU.

In the `energy0` stage the external forces are switched off. `--quick` runs a subset with shorter
runs (~30 s). `--smoke` runs 8 cases with `energy0` plus a 2000-step `nve`; it is meant for a fail-fast
check of a new (CUDA) build.

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

1. **AC with several chains (bug at f802bea4, fixed in anm_refactor 2dd5796a: `kdpg/AC` now passes).** `ANMInteraction::pair_interaction_bonded` computes
   `_computed_r = q->pos - p->pos` without the minimum image. oxDNA shifts every strand's centre of
   mass into the box when it reads a configuration (and periodically with `fix_diffusion`), so for a
   protein made of several chains the inter-chain springs see distances off by a box vector.
   `kdpg/AC`: step-0 energy 1720.74 per particle (original and numpy reference: 0.0821), the run
   blows up to NaN. Single-chain proteins and DNANM/DNACT (which use `min_image`) are not affected;
   `fix_diffusion = 0` hides it. The `test/DNANM/new_topology/kdpg_AC` reference value in oxDNA
   (2164.606502) was this wrong number (corrected in 52f59157).
2. **Protein torque from protein-nucleotide excluded volume (DNACT/RNACT).** The original adds the
   nucleotide's excluded-volume torque to the protein particle as well
   (`protein->torque -= protein->orientationT * -torquenuc`); the current code adds
   `(-r_hat * sigma/2) x F` about the protein centre, which is identically zero because the force is
   central (and uses `p->orientationT` instead of `protein->orientationT`, harmless only because the
   term is zero). The current behaviour conserves angular momentum, the original does not. Since
   bbbe1c00 the term is expressed in the amino acid's frame. It is still the near-zero central-force
   torque: the original's nucleotide torque was deliberately not reproduced. Effect: identical
   energies, and ANMT trajectories diverge at the first protein-nucleotide contact (`ubq_dna/DNACT`
   and `cage/DNACT` det are XFAIL; the deviation starts on the residue in contact and its backbone
   neighbours). Averages (stat) agree within error. The NVE check shows that the original is not
   energy conserving: over 20 time units the conserved energy drifts by +7.2e-3 per particle for
   `ubq_dna/DNACT` (original) vs -4.3e-5 (current), and by +2.3e-3 vs -1.6e-5 for `cage/DNACT`. The DNANM
   versions conserve energy in both codes.
3. **Debye-Hueckel in RNANM/RNACT.** The double counting of the DH term for RNA-RNA pairs was in the
   *pre-refactor current* oxDNA (`RNANMInteraction::pair_interaction_nonbonded` added `_debye_huckel`
   on top of `RNA2Interaction::pair_interaction_nonbonded`, which already includes it; removed in
   oxDNA commit 1f9f4e75). The original `anm-oxdna` RNANM/RNACT count it once. The `ubq_rna` cases
   measure it: the DH energy (E - E with DH off) is identical in both codes (ratio 1.000), and the
   total energies, trajectories and averages agree.
4. **Precision of the original ANMT.** See the det tolerances above (`sqrtf` in `normalize`).
5. **Parameter-file parsing.** At f802bea4 a malformed ANMT backbone line was silently skipped
   (`catch(...) { continue; }`); fixed in 8929388c (line-by-line parser).
6. **Input formats.** The original needs the protein neighbour lists in the topology; the current
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
  equilibrated.dat      start of the stat and nve runs (tools/equilibrate.py)
  stored.json           stored reference values ("old_code": "stored" cases)
  topology.top          new-style (5->3, type=...) topology
  anm.par / anmt.par    ANM and ANMT parameter files
  anm.top / mixed.top / hybrid.top   old-style topologies with neighbour lists
  anm_msf.dat           ANM-predicted fluctuations (protein-only systems)
tools/                  generators (see below)
scripts/                build scripts; cluster_cuda.sh = build + run everything on a GPU node
cuda_tests.md           instructions for the GPU runs
reference_results/      summaries of full CPU runs
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
# new-style topology versions (and the same for perturbed.dat / equilibrated.dat)
python tools/to_new_topology.py systems/ubq_dna/hybrid.top systems/ubq_dna/conf.dat systems/ubq_dna_nt --types DNA --forces systems/ubq_dna/forces.txt
python tools/to_new_topology.py systems/ubq_rna/hybrid.top systems/ubq_rna/conf.dat systems/ubq_rna_nt --types RNA --forces systems/ubq_rna/forces.txt
# ubiquitin + DNA/RNA hybrid duplex (na_init.* = oxDNA test/CUDA/NA/hybrid_overhang/init.*)
(cd tools && python make_hybrid_new.py ../systems/ubiquitin ../systems/ubq_hybrid/na_init.dat ../systems/ubq_hybrid/na_init.top ../systems/ubq_hybrid)
python tools/equilibrate.py <current oxDNA binary> ubq_hybrid
# na_*: init.dat/init.top of oxDNA test/CUDA/NA/<system> (as conf.dat/topology.top), perturbed with --sigma 0.01 --sigma-rot 0.05
```

`tools/models.py` is a copy of `ANMUtils/models.py` with three changes, all marked `PATCHED`:
`float(np.ravel(1 / slope)[0])` (numpy >= 2 refuses `float()` of a 1-element array), a headless
matplotlib backend, and normalised a3 vectors (the original writes `a3 = a1 x e_x` without
normalising it, while oxDNA normalises orientations on input, so the torsion reference
`d0 = a3_i . a3_j` in the par file did not correspond to the simulated configuration).

## Adding tests

- A new case: add an entry to `cases` in a system's `config.json`. Keys: `old`, `new` and `common`
  (input-file keys for each code / both; relative file names are resolved in the system directory),
  `old@<backend>`/`new@<backend>` (substitutions for one backend, e.g. `new@CUDA`), `old_code`,
  `reference` (par file and kb/kt for the numpy check), `quick`, `smoke`, `modes` (`original`,
  `same_code`), `stages_enabled[@mode]`, `cuda_only`, `xfail`
  (`{"<check>": reason, "<check>@CUDA_double": {"reason": ..., "match": "text the error must contain"},
  "*@CUDA": reason, "*@same_code": reason}`), `energy0_variants` / `stat_variants` (extra runs with
  modified inputs).
- A new system: a directory with the input files and a `config.json` (copy one). `protein_range` and
  `nucleic_range` select the particles for the RMSF; `cost` (ms per step) only orders the jobs.
- **Systems without an original-code reference** (oxNA, new interactions): use `"old_code": "stored"`,
  run once with `--update-stored`, and commit `stored.json`. Set `protein_range` to the particles
  whose RMSF should be compared (for the `na_*` systems, a duplex).

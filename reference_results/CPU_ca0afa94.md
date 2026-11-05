# ANM/oxNA regression results (CPU)

- mode: original vs current
- old/reference: `anm-oxdna (015a855) build, bin/oxDNA` (CPU)
- new/test: `oxDNA anm_refactor ca0afa94 build, bin/oxDNA` (CPU)
- wall time: 704 s

| case | check | status | details |
|---|---|---|---|
| cage/DNANM | inputs | PASS | 12720 springs, 12720 topology bonds, 0 springs not in topology |
| cage/DNANM | energy0 | PASS | default: ref -0.298704 test -0.298704 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.380868 test -0.380868 |
| cage/DNANM | det | PASS | 6 frames over 300 steps: max /dx/ 5.4e-04 (final 4.5e-04), max /dv,dL/ 4.1e-03, max /dE/ 3.0e-06 |
| cage/DNANM | stat | PASS | <E> old -0.50486+-0.00165 new -0.50405+-0.00140 (z +0.4); RMSF mean old 0.602 new 0.593 A, 0/605 residues /z/>4, <z^2> 1.20; RMSF(NA) mean old 3.503 new 3.743 A, 1/393 residues /z/>4, <z^2> 1.62 (not tested) |
| cage/DNANM | nve | PASS | conserved E drift over 20 time units: ref -1.62e-06 test -2.43e-05; fluct ref 1.7e-05 test 1.2e-05 (per particle) |
| cage/DNACT | inputs | PASS | 12720 springs, 12720 topology bonds, 0 springs not in topology |
| cage/DNACT | energy0 | PASS | default: ref -0.243505 test -0.243505 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.380868 test -0.380868 |
| cage/DNACT | det | XFAIL | 6 frames over 300 steps: max /dx/ 1.5e-02 (final 1.5e-02), max /dv,dL/ 2.5e-01, max /dE/ 1.6e-04; exceeds 0.002 from frame 1 (t = 100) [expected: by design: the original code also gives the amino acid the nucleotide's excluded-volume torque (protein->torque -= protein->orientationT * -torquenuc), which does not conserve angular momentum or energy (see the nve check of this case: the original drifts); the current code applies only (r x F) about the amino-acid centre in the amino acid's frame (bbbe1c00), which is ~0 for this central force, so ANMT trajectories diverge at the first protein-nucleotide contact] |
| cage/DNACT | stat | PASS | <E> old -0.43672+-0.00164 new -0.43496+-0.00132 (z +0.8); RMSF mean old 0.582 new 0.573 A, 2/605 residues /z/>4, <z^2> 1.51; RMSF(NA) mean old 3.494 new 3.804 A, 12/393 residues /z/>4, <z^2> 2.95 (not tested) |
| cage/DNACT | nve | PASS | conserved E drift over 20 time units: ref +2.33e-03 test -1.60e-05; fluct ref 3.6e-04 test 1.7e-05 (per particle) |
| cln025/AC | inputs | PASS | 29 springs, 29 topology bonds, 0 springs not in topology |
| cln025/AC | energy0 | PASS | default: ref 0.246245 test 0.246245 diff +0.00e+00; old-numpy +3.14e-07; new-numpy +3.14e-07 |
| cln025/AC | det | PASS | 10 frames over 5000 steps: max /dx/ 1.0e-14 (final 1.0e-14), max /dv,dL/ 2.5e-14, max /dE/ 0.0e+00 |
| cln025/AC | stat | PASS | <E> old 0.12725+-0.00407 new 0.12593+-0.00148 (z -0.3); RMSF mean old 0.721 new 0.692 A, 0/10 residues /z/>4, <z^2> 0.08; vs ANM prediction: r=0.62, sim/ANM=1.48 |
| cln025/AC | nve | PASS | conserved E drift over 20 time units: ref +1.51e-06 test +1.51e-06; fluct ref 1.6e-06 test 1.6e-06 (per particle) |
| cln025/DNACT | inputs | PASS | 29 springs, 29 topology bonds, 0 springs not in topology |
| cln025/DNACT | energy0 | PASS | default: ref 1.082646 test 1.082646 diff +0.00e+00; old-numpy -3.70e-08; new-numpy -3.70e-08 |
| cln025/DNACT | det | PASS | 10 frames over 5000 steps: max /dx/ 4.5e-08 (final 4.5e-08), max /dv,dL/ 1.4e-07, max /dE/ 0.0e+00 |
| cln025/DNACT | stat | PASS | <E> old 0.25267+-0.00413 new 0.26110+-0.00461 (z +1.4); RMSF mean old 0.395 new 0.402 A, 0/10 residues /z/>4, <z^2> 1.76; vs ANM prediction: r=0.94, sim/ANM=0.86 |
| cln025/DNACT | nve | PASS | conserved E drift over 20 time units: ref +1.92e-06 test +1.87e-06; fluct ref 2.7e-06 test 2.7e-06 (per particle) |
| cln025/RNACT | inputs | PASS | 29 springs, 29 topology bonds, 0 springs not in topology |
| cln025/RNACT | energy0 | PASS | default: ref 1.082646 test 1.082646 diff +0.00e+00; old-numpy -3.70e-08; new-numpy -3.70e-08 |
| cln025/RNACT | det | PASS | 10 frames over 5000 steps: max /dx/ 4.5e-08 (final 4.5e-08), max /dv,dL/ 1.4e-07, max /dE/ 0.0e+00 |
| cln025/RNACT | stat | PASS | <E> old 0.25267+-0.00413 new 0.26110+-0.00461 (z +1.4); RMSF mean old 0.395 new 0.402 A, 0/10 residues /z/>4, <z^2> 1.76; vs ANM prediction: r=0.94, sim/ANM=0.86 |
| cln025/RNACT | nve | PASS | conserved E drift over 20 time units: ref +1.92e-06 test +1.87e-06; fluct ref 2.7e-06 test 2.7e-06 (per particle) |
| kdpg/AC | inputs | PASS | 12716 springs, 12716 topology bonds, 0 springs not in topology |
| kdpg/AC | energy0 | PASS | default: ref 0.082090 test 0.082090 diff +0.00e+00; old-numpy -2.75e-07; new-numpy -2.75e-07 |
| kdpg/AC | det | PASS | 5 frames over 1000 steps: max /dx/ 2.0e-14 (final 2.0e-14), max /dv,dL/ 3.2e-13, max /dE/ 0.0e+00 |
| kdpg/AC | stat | PASS | <E> old 0.14884+-0.00111 new 0.14730+-0.00099 (z -1.0); RMSF mean old 0.567 new 0.560 A, 3/605 residues /z/>4, <z^2> 1.38 |
| kdpg/AC | nve | PASS | conserved E drift over 20 time units: ref +2.16e-04 test +2.15e-04; fluct ref 2.4e-05 test 2.4e-05 (per particle) |
| kdpg/ANMT | inputs | PASS | 12716 springs, 12716 topology bonds, 0 springs not in topology |
| kdpg/ANMT | energy0 | PASS | default: ref 1.209537 test 1.209537 diff +0.00e+00; old-numpy -2.69e-08; new-numpy -2.69e-08 |
| kdpg/ANMT | det | PASS | 5 frames over 1000 steps: max /dx/ 3.5e-07 (final 3.5e-07), max /dv,dL/ 1.4e-05, max /dE/ 0.0e+00 |
| kdpg/ANMT | stat | PASS | <E> old 0.31346+-0.00308 new 0.31360+-0.00286 (z +0.0); RMSF mean old 0.546 new 0.550 A, 4/605 residues /z/>4, <z^2> 1.47 |
| kdpg/ANMT | nve | PASS | conserved E drift over 20 time units: ref +4.13e-06 test -6.97e-05; fluct ref 2.6e-05 test 2.6e-05 (per particle) |
| kdpg/DNANM | inputs | PASS | 12716 springs, 12716 topology bonds, 0 springs not in topology |
| kdpg/DNANM | energy0 | PASS | default: ref 0.082090 test 0.082090 diff +0.00e+00; old-numpy -2.75e-07; new-numpy -2.75e-07 |
| kdpg/DNANM | det | PASS | 5 frames over 1000 steps: max /dx/ 2.0e-14 (final 2.0e-14), max /dv,dL/ 3.2e-13, max /dE/ 0.0e+00 |
| kdpg/DNANM | stat | PASS | <E> old 0.14884+-0.00111 new 0.14730+-0.00099 (z -1.0); RMSF mean old 0.567 new 0.560 A, 3/605 residues /z/>4, <z^2> 1.38 |
| kdpg/DNANM | nve | PASS | conserved E drift over 20 time units: ref +2.16e-04 test +2.15e-04; fluct ref 2.4e-05 test 2.4e-05 (per particle) |
| kdpg/ACT | inputs | PASS | 12716 springs, 12716 topology bonds, 0 springs not in topology |
| kdpg/ACT | energy0 | PASS | default: ref 1.209537 test 1.209537 diff +0.00e+00; old-numpy -2.69e-08; new-numpy -2.69e-08 |
| kdpg/ACT | det | PASS | 5 frames over 1000 steps: max /dx/ 3.5e-07 (final 3.5e-07), max /dv,dL/ 1.4e-05, max /dE/ 0.0e+00 |
| kdpg/ACT | stat | PASS | <E> old 0.31346+-0.00308 new 0.31360+-0.00286 (z +0.0); RMSF mean old 0.546 new 0.550 A, 4/605 residues /z/>4, <z^2> 1.47 |
| kdpg/ACT | nve | PASS | conserved E drift over 20 time units: ref +4.13e-06 test -6.97e-05; fluct ref 2.6e-05 test 2.6e-05 (per particle) |
| na_dna_duplex/NA | energy0 | PASS | default: ref -1.161387 test -1.161387 diff +0.00e+00 (stored ref) |
| na_dna_duplex/NA | stat | PASS | <E> old -1.27589+-0.01489 new -1.30450+-0.01084 (z -1.6); RMSF mean old 1.811 new 1.939 A, 0/22 residues /z/>4, <z^2> 0.83 |
| na_dna_duplex/NA | nve | PASS | conserved E drift over 20 time units: ref +1.30e-04 test +1.30e-04; fluct ref 5.8e-05 test 5.8e-05 (per particle) |
| na_hybrid_overhang/NA | energy0 | PASS | default: ref -1.295012 test -1.295012 diff +0.00e+00 (stored ref) |
| na_hybrid_overhang/NA | stat | PASS | <E> old -1.36194+-0.01634 new -1.35724+-0.01235 (z +0.2); RMSF mean old 2.499 new 2.988 A, 0/30 residues /z/>4, <z^2> 2.34 |
| na_hybrid_overhang/NA | nve | PASS | conserved E drift over 20 time units: ref +2.36e-05 test +2.36e-05; fluct ref 1.1e-05 test 1.1e-05 (per particle) |
| na_mixed/NA | energy0 | PASS | default: ref -1.186885 test -1.186885 diff +0.00e+00 (stored ref) |
| na_mixed/NA | stat | PASS | <E> old -1.25029+-0.00632 new -1.24846+-0.00317 (z +0.3); RMSF mean old 1.475 new 1.651 A, 0/22 residues /z/>4, <z^2> 0.78 |
| na_mixed/NA | nve | PASS | conserved E drift over 20 time units: ref -1.61e-05 test -1.61e-05; fluct ref 6.9e-06 test 6.9e-06 (per particle) |
| na_rna_duplex/NA | energy0 | PASS | default: ref -1.358674 test -1.358674 diff +0.00e+00 (stored ref) |
| na_rna_duplex/NA | stat | PASS | <E> old -1.46963+-0.00778 new -1.48338+-0.01177 (z -1.0); RMSF mean old 1.410 new 1.521 A, 0/24 residues /z/>4, <z^2> 1.30 |
| na_rna_duplex/NA | nve | PASS | conserved E drift over 20 time units: ref -1.50e-04 test -1.50e-04; fluct ref 3.8e-05 test 3.8e-05 (per particle) |
| ubiquitin/AC | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/AC | energy0 | PASS | default: ref 0.060555 test 0.060555 diff +0.00e+00; old-numpy +7.50e-08; new-numpy +7.50e-08 |
| ubiquitin/AC | det | PASS | 10 frames over 2000 steps: max /dx/ 1.0e-14 (final 1.0e-14), max /dv,dL/ 3.3e-14, max /dE/ 0.0e+00 |
| ubiquitin/AC | stat | PASS | <E> old 0.14589+-0.00039 new 0.14421+-0.00063 (z -2.3); RMSF mean old 0.637 new 0.627 A, 0/76 residues /z/>4, <z^2> 1.52; vs ANM prediction: r=0.98, sim/ANM=0.96 |
| ubiquitin/AC | nve | PASS | conserved E drift over 20 time units: ref -1.91e-04 test -1.91e-04; fluct ref 6.2e-05 test 6.2e-05 (per particle) |
| ubiquitin/ANMT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/ANMT | energy0 | PASS | default: ref 0.169886 test 0.169886 diff +0.00e+00; old-numpy -2.13e-07; new-numpy -2.13e-07 |
| ubiquitin/ANMT | det | PASS | 10 frames over 2000 steps: max /dx/ 1.3e-08 (final 1.3e-08), max /dv,dL/ 1.1e-07, max /dE/ 0.0e+00 |
| ubiquitin/ANMT | stat | PASS | <E> old 0.28389+-0.00576 new 0.28635+-0.00178 (z +0.4); RMSF mean old 0.553 new 0.557 A, 0/76 residues /z/>4, <z^2> 1.42; vs ANM prediction: r=0.99, sim/ANM=0.85 |
| ubiquitin/ANMT | nve | PASS | conserved E drift over 20 time units: ref -9.53e-05 test -9.55e-05; fluct ref 3.8e-05 test 3.8e-05 (per particle) |
| ubiquitin/DNACT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/DNACT | energy0 | PASS | default: ref 0.169886 test 0.169886 diff +0.00e+00; old-numpy -2.13e-07; new-numpy -2.13e-07 |
| ubiquitin/DNACT | det | PASS | 10 frames over 2000 steps: max /dx/ 1.3e-08 (final 1.3e-08), max /dv,dL/ 1.1e-07, max /dE/ 0.0e+00 |
| ubiquitin/DNACT | stat | PASS | <E> old 0.28322+-0.00441 new 0.28635+-0.00178 (z +0.7); RMSF mean old 0.553 new 0.557 A, 1/76 residues /z/>4, <z^2> 1.55; vs ANM prediction: r=0.99, sim/ANM=0.85 |
| ubiquitin/DNACT | nve | PASS | conserved E drift over 20 time units: ref -9.54e-05 test -9.55e-05; fluct ref 3.8e-05 test 3.8e-05 (per particle) |
| ubiquitin/DNANM | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/DNANM | energy0 | PASS | default: ref 0.060555 test 0.060555 diff +0.00e+00; old-numpy +7.50e-08; new-numpy +7.50e-08 |
| ubiquitin/DNANM | det | PASS | 10 frames over 2000 steps: max /dx/ 1.0e-14 (final 1.0e-14), max /dv,dL/ 3.3e-14, max /dE/ 0.0e+00 |
| ubiquitin/DNANM | stat | PASS | <E> old 0.14589+-0.00039 new 0.14421+-0.00063 (z -2.3); RMSF mean old 0.637 new 0.627 A, 0/76 residues /z/>4, <z^2> 1.52; vs ANM prediction: r=0.98, sim/ANM=0.96 |
| ubiquitin/DNANM | nve | PASS | conserved E drift over 20 time units: ref -1.91e-04 test -1.91e-04; fluct ref 6.2e-05 test 6.2e-05 (per particle) |
| ubiquitin/RNANM | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/RNANM | energy0 | PASS | default: ref 0.060555 test 0.060555 diff +0.00e+00; old-numpy +7.50e-08; new-numpy +7.50e-08 |
| ubiquitin/RNANM | det | PASS | 10 frames over 2000 steps: max /dx/ 1.0e-14 (final 1.0e-14), max /dv,dL/ 3.3e-14, max /dE/ 0.0e+00 |
| ubiquitin/RNANM | stat | PASS | <E> old 0.14589+-0.00039 new 0.14421+-0.00063 (z -2.3); RMSF mean old 0.637 new 0.627 A, 0/76 residues /z/>4, <z^2> 1.52; vs ANM prediction: r=0.98, sim/ANM=0.96 |
| ubiquitin/RNANM | nve | PASS | conserved E drift over 20 time units: ref -1.91e-04 test -1.91e-04; fluct ref 6.2e-05 test 6.2e-05 (per particle) |
| ubiquitin/RNACT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/RNACT | energy0 | PASS | default: ref 0.169886 test 0.169886 diff +0.00e+00; old-numpy -2.13e-07; new-numpy -2.13e-07 |
| ubiquitin/RNACT | det | PASS | 10 frames over 2000 steps: max /dx/ 1.3e-08 (final 1.3e-08), max /dv,dL/ 1.1e-07, max /dE/ 0.0e+00 |
| ubiquitin/RNACT | stat | PASS | <E> old 0.28322+-0.00441 new 0.28635+-0.00178 (z +0.7); RMSF mean old 0.553 new 0.557 A, 1/76 residues /z/>4, <z^2> 1.55; vs ANM prediction: r=0.99, sim/ANM=0.85 |
| ubiquitin/RNACT | nve | PASS | conserved E drift over 20 time units: ref -9.54e-05 test -9.55e-05; fluct ref 3.8e-05 test 3.8e-05 (per particle) |
| ubiquitin/ACT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/ACT | energy0 | PASS | default: ref 0.169886 test 0.169886 diff +0.00e+00; old-numpy -2.13e-07; new-numpy -2.13e-07 |
| ubiquitin/ACT | det | PASS | 10 frames over 2000 steps: max /dx/ 1.3e-08 (final 1.3e-08), max /dv,dL/ 1.1e-07, max /dE/ 0.0e+00 |
| ubiquitin/ACT | stat | PASS | <E> old 0.28389+-0.00576 new 0.28635+-0.00178 (z +0.4); RMSF mean old 0.553 new 0.557 A, 0/76 residues /z/>4, <z^2> 1.42; vs ANM prediction: r=0.99, sim/ANM=0.85 |
| ubiquitin/ACT | nve | PASS | conserved E drift over 20 time units: ref -9.53e-05 test -9.55e-05; fluct ref 3.8e-05 test 3.8e-05 (per particle) |
| ubq_dna/DNANM | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubq_dna/DNANM | energy0 | PASS | default: ref -0.241033 test -0.241033 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.287055 test -0.287055 |
| ubq_dna/DNANM | det | PASS | 10 frames over 500 steps: max /dx/ 1.3e-05 (final 5.5e-06), max /dv,dL/ 3.4e-04, max /dE/ 1.0e-06 |
| ubq_dna/DNANM | stat | PASS | <E> old -0.23331+-0.00208 new -0.23849+-0.00265 (z -1.5); RMSF mean old 0.635 new 0.620 A, 0/76 residues /z/>4, <z^2> 1.07; RMSF(NA) mean old 1.383 new 1.296 A, 0/24 residues /z/>4, <z^2> 1.69 (not tested) |
| ubq_dna/DNANM | nve | PASS | conserved E drift over 20 time units: ref +8.52e-06 test -1.30e-04; fluct ref 5.7e-05 test 4.7e-05 (per particle) |
| ubq_dna/DNACT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubq_dna/DNACT | energy0 | PASS | default: ref -0.157942 test -0.157942 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.287056 test -0.287056 |
| ubq_dna/DNACT | det | XFAIL | 10 frames over 500 steps: max /dx/ 3.9e-03 (final 1.8e-03), max /dv,dL/ 1.3e-01, max /dE/ 1.3e-03; exceeds 0.0002 from frame 0 (t = 50) [expected: by design: the original code also gives the amino acid the nucleotide's excluded-volume torque (protein->torque -= protein->orientationT * -torquenuc), which does not conserve angular momentum or energy (see the nve check of this case: the original drifts); the current code applies only (r x F) about the amino-acid centre in the amino acid's frame (bbbe1c00), which is ~0 for this central force, so ANMT trajectories diverge at the first protein-nucleotide contact] |
| ubq_dna/DNACT | stat | PASS | <E> old -0.12657+-0.00320 new -0.12865+-0.00293 (z -0.5); RMSF mean old 0.553 new 0.559 A, 0/76 residues /z/>4, <z^2> 1.03; RMSF(NA) mean old 1.348 new 1.315 A, 0/24 residues /z/>4, <z^2> 0.35 (not tested) |
| ubq_dna/DNACT | nve | PASS | conserved E drift over 20 time units: ref +7.19e-03 test -4.30e-05; fluct ref 1.2e-03 test 3.7e-05 (per particle) |
| ubq_dna_nt/NANM | energy0 | PASS | default: ref -0.241033 test -0.241033 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.287055 test -0.287055 |
| ubq_dna_nt/NANM | det | PASS | 10 frames over 1000 steps: max /dx/ 0.0e+00 (final 0.0e+00), max /dv,dL/ 0.0e+00, max /dE/ 0.0e+00 (bit-identical output) |
| ubq_dna_nt/NACT | energy0 | PASS | default: ref -0.157942 test -0.157942 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.287056 test -0.287056 |
| ubq_dna_nt/NACT | det | PASS | 10 frames over 1000 steps: max /dx/ 0.0e+00 (final 0.0e+00), max /dv,dL/ 0.0e+00, max /dE/ 0.0e+00 (bit-identical output) |
| ubq_hybrid/NANM | energy0 | PASS | default: ref -0.256789 test -0.256789 diff +0.00e+00 (stored ref); non-protein terms (E - E_numpy_protein): ref -0.300206 test -0.300206 |
| ubq_hybrid/NANM | stat | PASS | <E> old -0.29363+-0.00241 new -0.28764+-0.00252 (z +1.7); RMSF mean old 0.631 new 0.639 A, 0/76 residues /z/>4, <z^2> 1.04; RMSF(NA) mean old 2.194 new 2.240 A, 0/30 residues /z/>4, <z^2> 0.33 (not tested) |
| ubq_hybrid/NANM | nve | PASS | conserved E drift over 20 time units: ref +4.17e-04 test +4.17e-04; fluct ref 3.2e-05 test 3.2e-05 (per particle) |
| ubq_hybrid/NACT | energy0 | PASS | default: ref -0.178401 test -0.178401 diff +0.00e+00 (stored ref); non-protein terms (E - E_numpy_protein): ref -0.300206 test -0.300206 |
| ubq_hybrid/NACT | stat | PASS | <E> old -0.18730+-0.00260 new -0.19052+-0.00406 (z -0.7); RMSF mean old 0.553 new 0.552 A, 0/76 residues /z/>4, <z^2> 1.03; RMSF(NA) mean old 2.522 new 2.428 A, 0/30 residues /z/>4, <z^2> 0.20 (not tested) |
| ubq_hybrid/NACT | nve | PASS | conserved E drift over 20 time units: ref -2.59e-06 test -2.59e-06; fluct ref 3.0e-05 test 3.0e-05 (per particle) |
| ubq_rna/RNANM | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubq_rna/RNANM | energy0 | PASS | default: ref -0.274667 test -0.274667 diff +0.00e+00; _dh_off: ref -0.275158 test -0.275158 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.320689 test -0.320689; DH term ref 0.000491 test 0.000491 ratio 1.000 |
| ubq_rna/RNANM | det | PASS | 10 frames over 500 steps: max /dx/ 2.5e-05 (final 1.1e-05), max /dv,dL/ 1.6e-04, max /dE/ 1.0e-06 |
| ubq_rna/RNANM | stat_dh_off | PASS | <E> old -0.24810+-0.00361 new -0.25367+-0.00148 (z -1.4); RMSF mean old 0.622 new 0.620 A, 0/76 residues /z/>4, <z^2> 0.83; RMSF(NA) mean old 1.663 new 1.587 A, 0/24 residues /z/>4, <z^2> 0.54 (not tested) |
| ubq_rna/RNANM | stat | PASS | <E> old -0.25260+-0.00215 new -0.25173+-0.00091 (z +0.4); RMSF mean old 0.621 new 0.622 A, 0/76 residues /z/>4, <z^2> 1.02; RMSF(NA) mean old 1.597 new 1.501 A, 0/24 residues /z/>4, <z^2> 1.30 (not tested) |
| ubq_rna/RNANM | nve | PASS | conserved E drift over 20 time units: ref -9.41e-05 test -1.48e-04; fluct ref 2.0e-05 test 2.1e-05 (per particle) |
| ubq_rna/RNACT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubq_rna/RNACT | energy0 | PASS | default: ref -0.191575 test -0.191575 diff +0.00e+00; _dh_off: ref -0.192066 test -0.192066 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.320689 test -0.320689; DH term ref 0.000491 test 0.000491 ratio 1.000 |
| ubq_rna/RNACT | det | PASS | 10 frames over 500 steps: max /dx/ 1.3e-05 (final 1.3e-05), max /dv,dL/ 2.4e-04, max /dE/ 1.0e-06 |
| ubq_rna/RNACT | stat_dh_off | PASS | <E> old -0.14162+-0.00279 new -0.14417+-0.00213 (z -0.7); RMSF mean old 0.552 new 0.561 A, 0/76 residues /z/>4, <z^2> 1.20; RMSF(NA) mean old 1.582 new 1.509 A, 0/24 residues /z/>4, <z^2> 0.57 (not tested) |
| ubq_rna/RNACT | stat | PASS | <E> old -0.14127+-0.00419 new -0.14271+-0.00125 (z -0.3); RMSF mean old 0.552 new 0.551 A, 0/76 residues /z/>4, <z^2> 0.96; RMSF(NA) mean old 1.506 new 1.524 A, 0/24 residues /z/>4, <z^2> 0.14 (not tested) |
| ubq_rna/RNACT | nve | PASS | conserved E drift over 20 time units: ref -1.07e-03 test -1.19e-04; fluct ref 2.8e-04 test 3.3e-05 (per particle) |
| ubq_rna_nt/NANM | energy0 | PASS | default: ref -0.274667 test -0.274667 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.320689 test -0.320689 |
| ubq_rna_nt/NANM | det | PASS | 10 frames over 1000 steps: max /dx/ 0.0e+00 (final 0.0e+00), max /dv,dL/ 0.0e+00, max /dE/ 0.0e+00 (bit-identical output) |
| ubq_rna_nt/NACT | energy0 | PASS | default: ref -0.191575 test -0.191575 diff +0.00e+00; non-protein terms (E - E_numpy_protein): ref -0.320689 test -0.320689 |
| ubq_rna_nt/NACT | det | PASS | 10 frames over 1000 steps: max /dx/ 0.0e+00 (final 0.0e+00), max /dv,dL/ 0.0e+00, max /dE/ 0.0e+00 (bit-identical output) |

totals: PASS 126, XFAIL 2

**0 unexpected failure(s)**

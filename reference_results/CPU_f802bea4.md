# ANM regression results (CPU)

- old: `anm-oxdna (015a855) build, bin/oxDNA` (CPU)
- new: `oxDNA f802bea4 build, bin/oxDNA` (CPU)
- wall time: 828 s

| case | check | status | details |
|---|---|---|---|
| cage/DNANM | inputs | PASS | 12720 springs, 12720 topology bonds, 0 springs not in topology |
| cage/DNANM | energy0 | PASS | default: old -0.298704 new -0.298704 diff +0.00e+00; non-protein terms (E - E_ref_protein): old -0.380868 new -0.380868 |
| cage/DNANM | det | PASS | 6 frames over 300 steps: max /dx/ 5.4e-04 (final 4.5e-04), max /dv,dL/ 4.1e-03, max /dE/ 3.0e-06 |
| cage/DNANM | stat | PASS | <E> old -0.50486+-0.00165 new -0.50422+-0.00134 (z +0.3); RMSF mean old 0.602 new 0.593 A, 1/605 residues /z/>4, <z^2> 1.24; RMSF(NA) mean old 3.503 new 3.741 A, 1/393 residues /z/>4, <z^2> 1.60 (not tested) |
| cage/DNACT | inputs | PASS | 12720 springs, 12720 topology bonds, 0 springs not in topology |
| cage/DNACT | energy0 | PASS | default: old -0.243505 new -0.243505 diff +0.00e+00; non-protein terms (E - E_ref_protein): old -0.380868 new -0.380868 |
| cage/DNACT | det | XFAIL | 6 frames over 300 steps: max /dx/ 1.5e-02 (final 1.5e-02), max /dv,dL/ 2.5e-01, max /dE/ 1.6e-04; exceeds 0.002 from frame 1 (t = 100) [expected: the original code adds the nucleotide's excluded-volume torque to the protein particle (protein->torque -= protein->orientationT * -torquenuc); the current code adds (r x F) about the protein centre, which is zero for this central force, so ANMT trajectories diverge as soon as a protein-nucleotide contact occurs] |
| cage/DNACT | stat | PASS | <E> old -0.43672+-0.00164 new -0.43492+-0.00134 (z +0.9); RMSF mean old 0.582 new 0.573 A, 1/605 residues /z/>4, <z^2> 1.44; RMSF(NA) mean old 3.494 new 3.800 A, 12/393 residues /z/>4, <z^2> 2.94 (not tested) |
| cln025/AC | inputs | PASS | 29 springs, 29 topology bonds, 0 springs not in topology |
| cln025/AC | energy0 | PASS | default: old 0.246245 new 0.246245 diff +0.00e+00; old-ref +3.14e-07; new-ref +3.14e-07 |
| cln025/AC | det | PASS | 10 frames over 5000 steps: max /dx/ 6.7e-14 (final 6.7e-14), max /dv,dL/ 7.6e-14, max /dE/ 0.0e+00 |
| cln025/AC | stat | PASS | <E> old 0.12725+-0.00407 new 0.12698+-0.00270 (z -0.1); RMSF mean old 0.721 new 0.657 A, 0/10 residues /z/>4, <z^2> 0.39; vs ANM prediction: r=0.68, sim/ANM=1.41 |
| cln025/DNACT | inputs | PASS | 29 springs, 29 topology bonds, 0 springs not in topology |
| cln025/DNACT | energy0 | PASS | default: old 1.082646 new 1.082646 diff +0.00e+00; old-ref -3.70e-08; new-ref -3.70e-08 |
| cln025/DNACT | det | PASS | 10 frames over 5000 steps: max /dx/ 4.5e-08 (final 4.5e-08), max /dv,dL/ 1.4e-07, max /dE/ 0.0e+00 |
| cln025/DNACT | stat | PASS | <E> old 0.25267+-0.00413 new 0.25816+-0.00303 (z +1.1); RMSF mean old 0.395 new 0.401 A, 0/10 residues /z/>4, <z^2> 0.88; vs ANM prediction: r=0.94, sim/ANM=0.86 |
| cln025/RNACT | inputs | PASS | 29 springs, 29 topology bonds, 0 springs not in topology |
| cln025/RNACT | energy0 | PASS | default: old 1.082646 new 1.082646 diff +0.00e+00; old-ref -3.70e-08; new-ref -3.70e-08 |
| cln025/RNACT | det | PASS | 10 frames over 5000 steps: max /dx/ 4.5e-08 (final 4.5e-08), max /dv,dL/ 1.4e-07, max /dE/ 0.0e+00 |
| cln025/RNACT | stat | PASS | <E> old 0.25267+-0.00413 new 0.25816+-0.00303 (z +1.1); RMSF mean old 0.395 new 0.401 A, 0/10 residues /z/>4, <z^2> 0.88; vs ANM prediction: r=0.94, sim/ANM=0.86 |
| kdpg/AC | inputs | PASS | 12716 springs, 12716 topology bonds, 0 springs not in topology |
| kdpg/AC | energy0 | FAIL | default: old 0.082090 new 1720.743625 diff +1.72e+03; old-ref -2.75e-07; new-ref +1.72e+03 |
| kdpg/AC | det | FAIL | 5 frames over 1000 steps: max /dx/ 1.1e+01 (final 1.1e+01), max /dv,dL/ 1.0e+15, max /dE/ 4.1e+27; exceeds 1e-06 from frame 0 (t = 200) |
| kdpg/AC | stat | FAIL | rep0: non-finite or missing energies |
| kdpg/ANMT | inputs | PASS | 12716 springs, 12716 topology bonds, 0 springs not in topology |
| kdpg/ANMT | energy0 | PASS | default: old 1.209537 new 1.209537 diff +0.00e+00; old-ref -2.69e-08; new-ref -2.69e-08 |
| kdpg/ANMT | det | PASS | 5 frames over 1000 steps: max /dx/ 3.5e-07 (final 3.5e-07), max /dv,dL/ 1.4e-05, max /dE/ 0.0e+00 |
| kdpg/ANMT | stat | PASS | <E> old 0.31346+-0.00308 new 0.31360+-0.00254 (z +0.0); RMSF mean old 0.546 new 0.549 A, 3/605 residues /z/>4, <z^2> 1.45 |
| kdpg/DNANM | inputs | PASS | 12716 springs, 12716 topology bonds, 0 springs not in topology |
| kdpg/DNANM | energy0 | PASS | default: old 0.082090 new 0.082090 diff +0.00e+00; old-ref -2.75e-07; new-ref -2.75e-07 |
| kdpg/DNANM | det | PASS | 5 frames over 1000 steps: max /dx/ 5.0e-13 (final 3.0e-13), max /dv,dL/ 4.5e-11, max /dE/ 0.0e+00 |
| kdpg/DNANM | stat | PASS | <E> old 0.14884+-0.00111 new 0.14712+-0.00082 (z -1.2); RMSF mean old 0.567 new 0.558 A, 4/605 residues /z/>4, <z^2> 1.48 |
| ubiquitin/AC | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/AC | energy0 | PASS | default: old 0.060555 new 0.060555 diff +0.00e+00; old-ref +7.50e-08; new-ref +7.50e-08 |
| ubiquitin/AC | det | PASS | 10 frames over 2000 steps: max /dx/ 1.0e-14 (final 1.0e-14), max /dv,dL/ 7.3e-14, max /dE/ 0.0e+00 |
| ubiquitin/AC | stat | PASS | <E> old 0.14589+-0.00039 new 0.14376+-0.00195 (z -1.1); RMSF mean old 0.637 new 0.620 A, 2/76 residues /z/>4, <z^2> 1.98; vs ANM prediction: r=0.98, sim/ANM=0.95 |
| ubiquitin/ANMT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/ANMT | energy0 | PASS | default: old 0.169886 new 0.169886 diff +0.00e+00; old-ref -2.13e-07; new-ref -2.13e-07 |
| ubiquitin/ANMT | det | PASS | 10 frames over 2000 steps: max /dx/ 1.3e-08 (final 1.3e-08), max /dv,dL/ 1.1e-07, max /dE/ 0.0e+00 |
| ubiquitin/ANMT | stat | PASS | <E> old 0.28389+-0.00576 new 0.28542+-0.00287 (z +0.2); RMSF mean old 0.553 new 0.547 A, 1/76 residues /z/>4, <z^2> 1.21; vs ANM prediction: r=0.99, sim/ANM=0.84 |
| ubiquitin/DNACT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/DNACT | energy0 | PASS | default: old 0.169886 new 0.169886 diff +0.00e+00; old-ref -2.13e-07; new-ref -2.13e-07 |
| ubiquitin/DNACT | det | PASS | 10 frames over 2000 steps: max /dx/ 1.3e-08 (final 1.3e-08), max /dv,dL/ 1.1e-07, max /dE/ 0.0e+00 |
| ubiquitin/DNACT | stat | PASS | <E> old 0.28322+-0.00441 new 0.28542+-0.00287 (z +0.4); RMSF mean old 0.553 new 0.547 A, 0/76 residues /z/>4, <z^2> 1.39; vs ANM prediction: r=0.99, sim/ANM=0.84 |
| ubiquitin/DNANM | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/DNANM | energy0 | PASS | default: old 0.060555 new 0.060555 diff +0.00e+00; old-ref +7.50e-08; new-ref +7.50e-08 |
| ubiquitin/DNANM | det | PASS | 10 frames over 2000 steps: max /dx/ 1.0e-14 (final 1.0e-14), max /dv,dL/ 7.3e-14, max /dE/ 0.0e+00 |
| ubiquitin/DNANM | stat | PASS | <E> old 0.14589+-0.00039 new 0.14376+-0.00195 (z -1.1); RMSF mean old 0.637 new 0.620 A, 2/76 residues /z/>4, <z^2> 1.98; vs ANM prediction: r=0.98, sim/ANM=0.95 |
| ubiquitin/RNANM | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/RNANM | energy0 | PASS | default: old 0.060555 new 0.060555 diff +0.00e+00; old-ref +7.50e-08; new-ref +7.50e-08 |
| ubiquitin/RNANM | det | PASS | 10 frames over 2000 steps: max /dx/ 1.0e-14 (final 1.0e-14), max /dv,dL/ 7.3e-14, max /dE/ 0.0e+00 |
| ubiquitin/RNANM | stat | PASS | <E> old 0.14589+-0.00039 new 0.14376+-0.00195 (z -1.1); RMSF mean old 0.637 new 0.620 A, 2/76 residues /z/>4, <z^2> 1.98; vs ANM prediction: r=0.98, sim/ANM=0.95 |
| ubiquitin/RNACT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubiquitin/RNACT | energy0 | PASS | default: old 0.169886 new 0.169886 diff +0.00e+00; old-ref -2.13e-07; new-ref -2.13e-07 |
| ubiquitin/RNACT | det | PASS | 10 frames over 2000 steps: max /dx/ 1.3e-08 (final 1.3e-08), max /dv,dL/ 1.1e-07, max /dE/ 0.0e+00 |
| ubiquitin/RNACT | stat | PASS | <E> old 0.28322+-0.00441 new 0.28542+-0.00287 (z +0.4); RMSF mean old 0.553 new 0.547 A, 0/76 residues /z/>4, <z^2> 1.39; vs ANM prediction: r=0.99, sim/ANM=0.84 |
| ubq_dna/DNANM | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubq_dna/DNANM | energy0 | PASS | default: old -0.241033 new -0.241033 diff +0.00e+00; non-protein terms (E - E_ref_protein): old -0.287055 new -0.287055 |
| ubq_dna/DNANM | det | PASS | 10 frames over 500 steps: max /dx/ 1.3e-05 (final 5.5e-06), max /dv,dL/ 3.4e-04, max /dE/ 1.0e-06 |
| ubq_dna/DNANM | stat | PASS | <E> old -0.23331+-0.00208 new -0.23720+-0.00240 (z -1.2); RMSF mean old 0.635 new 0.626 A, 2/76 residues /z/>4, <z^2> 1.49; RMSF(NA) mean old 1.383 new 1.307 A, 0/24 residues /z/>4, <z^2> 2.03 (not tested) |
| ubq_dna/DNACT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubq_dna/DNACT | energy0 | PASS | default: old -0.157942 new -0.157942 diff +0.00e+00; non-protein terms (E - E_ref_protein): old -0.287056 new -0.287056 |
| ubq_dna/DNACT | det | XFAIL | 10 frames over 500 steps: max /dx/ 3.9e-03 (final 1.8e-03), max /dv,dL/ 1.3e-01, max /dE/ 1.3e-03; exceeds 0.0002 from frame 0 (t = 50) [expected: the original code adds the nucleotide's excluded-volume torque to the protein particle (protein->torque -= protein->orientationT * -torquenuc); the current code adds (r x F) about the protein centre, which is zero for this central force, so ANMT trajectories diverge as soon as a protein-nucleotide contact occurs] |
| ubq_dna/DNACT | stat | PASS | <E> old -0.12657+-0.00320 new -0.13287+-0.00245 (z -1.6); RMSF mean old 0.553 new 0.560 A, 0/76 residues /z/>4, <z^2> 1.14; RMSF(NA) mean old 1.348 new 1.280 A, 0/24 residues /z/>4, <z^2> 1.24 (not tested) |
| ubq_rna/RNANM | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubq_rna/RNANM | energy0 | PASS | default: old -0.274667 new -0.274667 diff +0.00e+00; _dh_off: old -0.275158 new -0.275158 diff +0.00e+00; non-protein terms (E - E_ref_protein): old -0.320689 new -0.320689; DH term old 0.000491 new 0.000491 ratio 1.000 |
| ubq_rna/RNANM | det | PASS | 10 frames over 500 steps: max /dx/ 2.5e-05 (final 1.1e-05), max /dv,dL/ 1.6e-04, max /dE/ 1.0e-06 |
| ubq_rna/RNANM | stat_dh_off | PASS | <E> old -0.24810+-0.00361 new -0.25288+-0.00233 (z -1.1); RMSF mean old 0.622 new 0.617 A, 0/76 residues /z/>4, <z^2> 1.14; RMSF(NA) mean old 1.663 new 1.529 A, 0/24 residues /z/>4, <z^2> 1.43 (not tested) |
| ubq_rna/RNANM | stat | PASS | <E> old -0.25260+-0.00215 new -0.25067+-0.00242 (z +0.6); RMSF mean old 0.621 new 0.625 A, 0/76 residues /z/>4, <z^2> 1.27; RMSF(NA) mean old 1.597 new 1.542 A, 0/24 residues /z/>4, <z^2> 1.06 (not tested) |
| ubq_rna/RNACT | inputs | PASS | 1033 springs, 1033 topology bonds, 0 springs not in topology |
| ubq_rna/RNACT | energy0 | PASS | default: old -0.191575 new -0.191575 diff +0.00e+00; _dh_off: old -0.192066 new -0.192066 diff +0.00e+00; non-protein terms (E - E_ref_protein): old -0.320689 new -0.320689; DH term old 0.000491 new 0.000491 ratio 1.000 |
| ubq_rna/RNACT | det | PASS | 10 frames over 500 steps: max /dx/ 1.3e-05 (final 1.3e-05), max /dv,dL/ 2.4e-04, max /dE/ 1.0e-06 |
| ubq_rna/RNACT | stat_dh_off | PASS | <E> old -0.14162+-0.00279 new -0.14250+-0.00278 (z -0.2); RMSF mean old 0.552 new 0.562 A, 1/76 residues /z/>4, <z^2> 1.76; RMSF(NA) mean old 1.582 new 1.492 A, 0/24 residues /z/>4, <z^2> 0.76 (not tested) |
| ubq_rna/RNACT | stat | PASS | <E> old -0.14127+-0.00419 new -0.14509+-0.00198 (z -0.8); RMSF mean old 0.552 new 0.555 A, 0/76 residues /z/>4, <z^2> 1.02; RMSF(NA) mean old 1.506 new 1.430 A, 0/24 residues /z/>4, <z^2> 1.20 (not tested) |

**3 unexpected failure(s)**

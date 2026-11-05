#!/usr/bin/env python3
"""Write systems/<name>/equilibrated.dat: the last configuration of an MD run of the ANMT variant of the
system (orientations matter only for ANMT), started from conf.dat. The 'stat' runs of both codes start
from this file so that the averages are not dominated by the thermalisation transient.

usage: tools/equilibrate.py OXDNA_BINARY SYSTEM [SYSTEM...]
"""
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import run_tests  # noqa: E402

SETTINGS = {  # system: (case providing the inputs, steps)
    'ubiquitin': ('DNACT', 100000),
    'cln025': ('DNACT', 200000),
    'kdpg': ('ANMT', 30000),
    'ubq_dna': ('DNACT', 100000),
    'ubq_rna': ('RNACT', 100000),
    'cage': ('DNACT', 20000),
    'ubq_hybrid': ('NACT', 100000),
}


def main():
    binary = os.path.abspath(sys.argv[1])
    for name in sys.argv[2:]:
        case_name, steps = SETTINGS[name]
        sysdir = os.path.join(ROOT, 'systems', name)
        cfg = json.load(open(os.path.join(sysdir, 'config.json')))
        case = dict(cfg['cases'][case_name], system=name, name=f'{name}/{case_name}', dir=sysdir, sys=cfg)
        with tempfile.TemporaryDirectory() as tmp:
            run_tests.make_input(case, ['new'], 'CPU', {'steps': steps, 'seed': 4242, 'conf_file': 'conf.dat',
                                                      'print_energy_every': steps // 10,
                                                      'print_conf_interval': steps}, tmp)
            subprocess.run([binary, 'input'], cwd=tmp, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            with open(os.path.join(tmp, 'last_conf.dat')) as f:
                lines = f.readlines()
            lines[0] = 't = 0\n'
            with open(os.path.join(sysdir, 'equilibrated.dat'), 'w') as f:
                f.writelines(lines)
            print(name, open(os.path.join(tmp, 'energy.dat')).read().splitlines()[-1])


if __name__ == '__main__':
    main()

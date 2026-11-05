#!/usr/bin/env python3
"""Convert an old-style (3'->5') hybrid protein/nucleic-acid system to a new-style (5'->3') topology.

    to_new_topology.py OLD_TOP CONF OUT_DIR --types DNA [RNA ...] [--forces forces.txt]

Writes OUT_DIR/topology.top, conf.dat (reordered: every nucleic-acid strand is reversed) and, if given,
forces.txt with the particle indices remapped. Amino-acid strands keep their order (N to C terminus in
both formats), so parameter files are unchanged as long as the protein comes first. --types gives the
type (DNA or RNA) of each nucleic-acid strand in order of appearance (one value applies to all).
Peptide strands are written with type=peptide, so the result is valid for NANM/NACT as well as for
DNANM/DNACT/RNANM/RNACT.
"""
import argparse
import os
import re


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('top')
    ap.add_argument('conf')
    ap.add_argument('outdir')
    ap.add_argument('--types', nargs='+', required=True)
    ap.add_argument('--forces')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    with open(a.top) as f:
        header = f.readline().split()
        parts = [l.split() for l in f if l.strip()]
    N = int(header[0])
    # group particles by strand, keeping file order
    strands, order = {}, []
    for i, p in enumerate(parts):
        sid = int(p[0])
        if sid not in strands:
            strands[sid] = []
            order.append(sid)
        strands[sid].append((i, p[1]))
    lines, new_index = [], []
    nuc_k = 0
    for sid in order:
        members = strands[sid]
        if sid < 0:
            seq = ''.join(c for _, c in members)
            lines.append(f'{seq} type=peptide')
            new_index += [i for i, _ in members]
        else:
            t = a.types[min(nuc_k, len(a.types) - 1)].upper()
            nuc_k += 1
            members = members[::-1]  # 5' -> 3'
            seq = ''.join(c for _, c in members)
            if t == 'RNA':
                seq = seq.replace('T', 'U')
            lines.append(f'{seq} type={t}')
            new_index += [i for i, _ in members]
    assert len(new_index) == N
    old_to_new = {old: new for new, old in enumerate(new_index)}

    with open(os.path.join(a.outdir, 'topology.top'), 'w') as f:
        f.write(f'{N} {len(order)} 5->3\n')
        f.write('\n'.join(lines) + '\n')

    with open(a.conf) as f:
        clines = f.readlines()
    head, body = clines[:3], [l for l in clines[3:] if l.strip()]
    with open(os.path.join(a.outdir, 'conf.dat'), 'w') as f:
        f.writelines(head)
        for old in new_index:
            f.write(body[old])

    if a.forces:
        with open(a.forces) as f:
            txt = f.read()
        txt = re.sub(r'^(\s*(?:particle|ref_particle)\s*=\s*)(\d+)',
                     lambda m: m.group(1) + str(old_to_new[int(m.group(2))]), txt, flags=re.M)
        with open(os.path.join(a.outdir, 'forces.txt'), 'w') as f:
            f.write(txt)


if __name__ == '__main__':
    main()

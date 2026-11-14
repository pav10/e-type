#!/usr/bin/env python3
"""Group saved survivors (same field, same brick) by O_K-isometry of the Leech structure.
usage: group_survivors.py N t"""
import sys, os, re, json, random, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.represent import herm_isom, herm_aut
from etype.certify import ok_basis
N, t = int(sys.argv[1]), int(sys.argv[2]); K = (N, t); d = disc(K)
leech = load("leech_d%d.json" % d); random.seed(0); random.shuffle(leech)
bricks = json.load(open(os.path.join(RES, "bricks_d%d_m4.json" % d)))
surv = collections.defaultdict(list)
for f in sorted(os.listdir(RES)):
    m = re.match(r"SURVIVOR_rank96_d%d_\d+_\d+_\w+?_(\d+)_(\d+)\.json$" % d, f)
    if m:
        surv[int(m.group(1))].append(int(m.group(2)))
groups = []
for b, ls in surv.items():
    A, JA = pari(bricks[b]["G"]), pari(bricks[b]["J"])
    try:
        ok_basis(A, JA, K); free = True
    except RuntimeError:
        free = False
    reps = []
    for l in ls:
        P, JP = leech[l]["G"], leech[l]["J"]
        o = herm_aut(P, JP, K)[0]
        print("  brick #%d Leech#%d |Aut_h| = %d" % (b, l, o), flush=True)
        for g in reps:
            if g["aut"] == o and herm_isom(P, JP, leech[g["rep"]]["G"], leech[g["rep"]]["J"], K):
                g["members"].append(l); break
        else:
            reps.append({"brick": b, "rep": l, "aut": o, "members": [l], "brick_free": free})
    groups += reps
for g in groups:
    print("brick #%d (free: %s)  Leech rep #%d  |Aut_h(P)| = %d  members %s" % (g["brick"], g["brick_free"], g["rep"], g["aut"], g["members"]), flush=True)
json.dump(groups, open(os.path.join(RES, "survivor_groups_d%d.json" % d), "w"))

#!/usr/bin/env python3
"""Add shapes for non-free sections to a brick-shape cache made before steinitz_gram existed.
usage: patch_nonfree.py CACHE.pkl N t"""
import sys, os, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.exact import steinitz_gram, _shape_worker
f = sys.argv[1]; K = (int(sys.argv[2]), int(sys.argv[3]))
shapes, nonfree, v = pickle.load(open(f, "rb"))
new_nf = []; added = 0
for item in nonfree:
    if len(item) == 4:
        new_nf.append(item); continue
    r, Gs, Js = item
    W, ideals = steinitz_gram(pari(Gs), pari(Js), K)
    sh = [(r,) + s for s in _shape_worker((W, K, 4.0, 0, 1))]
    shapes += sh; added += len(sh)
    new_nf.append((r, Gs, Js, ideals))
pickle.dump((shapes, new_nf, v), open(f, "wb"))
print("non-free sections: %d, shapes added: %d, total shapes: %d" % (len(new_nf), added, len(shapes)))

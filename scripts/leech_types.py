#!/usr/bin/env python3
"""Complete table of Hermitian inner products between minimal vectors of the
Leech structures leech_d19.json[idx] (shuffled order used by rank96_tensor.py).
usage: leech_types.py idx [idx ...]"""
import sys, os, time, random, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
K = (5, 1)
leech = load("leech_d19.json"); random.seed(0); random.shuffle(leech)
for idx in map(int, sys.argv[1:]):
    t = time.time()
    P, JP = leech[idx]["G"], leech[idx]["J"]
    V = short_vectors(P, 4)
    import numpy as np
    V = np.concatenate([V, -V], axis=1)                # both signs: all 196560
    tp = pair_types(P, JP, K, V)
    pickle.dump(tp, open(os.path.join(RES, "types22_d19_%d.pkl" % idx), "wb"))
    print("Leech#%d: %d types among minimal vectors (%.0fs)" % (idx, len(tp), time.time() - t), flush=True)

#!/usr/bin/env python3
"""Structured kill test on the d=19 rank-96 survivors:  rank2_test.py [normP] [sample]"""
import sys, os, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
K = (5, 1)
normP = int(sys.argv[1]) if len(sys.argv) > 1 else 4
sample = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
A, JA = {k: v for k, v in brick_walk(K, 4, steps=2500).items()}[(6, 20)]
leech = load("leech_d19.json"); random.seed(0); random.shuffle(leech)
rows = load("rank96_tensor_d19.json")
surv = [r["leech"] for r in rows if r["shortest_found"] >= 10]
print("survivors:", surv, flush=True)
for i in surv + [0, 1]:                      # two dead ones as controls
    t = time.time()
    kills, best = rank2_kill(A, JA, leech[i]["G"], leech[i]["J"], K, 8, normA=8, normP=normP, sample=sample)
    print("Leech#%d (m=%d): %d killing section pairs, best %s   (%.0fs)" % (i, leech[i]["m"], len(kills), best, time.time() - t), kills[:2], flush=True)

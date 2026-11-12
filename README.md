# e-type

Search for explicit lattices that are not of E-type (Kitaoka) at small rank:
extremal even unimodular of rank 96, extremal 3-modular of rank 84, and the
related extension of Kitaoka's rank bound.

Start with **[PLAN.md](PLAN.md)** — critical review of the earlier notes, the
constructions, validation, first results and the next steps.

```bash
./setup.sh
```

```bash
python3 scripts/validate.py
```

```bash
python3 scripts/rank96_tensor.py 5 1
```

Python 3 with cypari2 and fpylll, plus PARI/GP code in `gp/`. No Magma needed.

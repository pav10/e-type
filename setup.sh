#!/bin/sh
# One-time setup: Python deps + fplll BKZ strategies (not shipped in the fpylll wheel).
set -e
pip install cypari2 fpylll numpy
mkdir -p data
[ -s data/fplll_default_strategies.json ] || curl -sSL -o data/fplll_default_strategies.json \
  https://raw.githubusercontent.com/fplll/fplll/master/strategies/default.json
echo "ok"

#!/usr/bin/env python3
"""run_all.py — canonical, portable entry point.

Regenerates every figure from cache or from scratch (PLAN §9). Currently a
scaffold: the pipeline is wired once the validation ladder (PLAN §4) passes and
docs/model_spec.md lands.

Usage:
    python run_all.py            # run the full pipeline (E1..E6, gated)
    python run_all.py --list     # list stages
"""
from __future__ import annotations

import argparse
import importlib

# (name, module, note). E4 is a STATED NULL RESULT (docs/decisions.md): the route
# to chaos has only 2 successive period-doublings, so the Feigenbaum test is not
# run — reporting delta from 2 doublings would overclaim. It is skipped here.
STAGES = [
    ("e1", "experiments.e1_regime_map", None),
    ("e2", "experiments.e2_lyapunov", "chaos GATE"),
    ("e3", "experiments.e3_route", None),
    ("e4", None, "SKIPPED — stated null (only 2 doublings; no Feigenbaum cascade)"),
    ("e5", "experiments.e5_sensitivity", None),
    ("e6", "experiments.e6_multistability", None),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list stages and exit")
    args = parser.parse_args()

    if args.list:
        for name, mod, note in STAGES:
            print(f"{name}: {mod or '(skipped)'}" + (f"  [{note}]" if note else ""))
        return

    for name, mod, note in STAGES:
        header = f"=== {name.upper()}" + (f" ({note})" if note else "") + " ==="
        print(f"\n{header}")
        if mod is None:
            print(f"  {note}")
            continue
        importlib.import_module(mod).main()   # each stage is idempotent (cached)


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()

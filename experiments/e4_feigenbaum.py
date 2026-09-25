"""E4 — Feigenbaum test (ONLY if >=3 successive period-doublings resolve).

See PLAN §5-E4. Locate doublings by Newton shooting + continuation, detecting
Floquet multiplier = -1 (not by eye). Hypothesis to TEST: delta_n -> 4.669.
Report delta_n +/- error with the limited-ratio caveat.
"""
from __future__ import annotations


def main(config_path: str = "config/e4.yaml") -> None:
    raise NotImplementedError("Scaffold only; see PLAN §5-E4.")


if __name__ == "__main__":
    main()

from __future__ import annotations


def parse_xtb_energy(text: str) -> float | None:
    for line in text.splitlines():
        if "TOTAL ENERGY" in line.upper():
            parts = line.replace("=", " ").split()
            for part in reversed(parts):
                try:
                    return float(part)
                except ValueError:
                    continue
    return None


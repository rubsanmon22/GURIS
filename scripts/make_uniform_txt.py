#!/usr/bin/env python3
import numpy as np
from pathlib import Path
import re
import sys


def strip_comment(line):
    """Remove Fortran namelist comments starting with !."""
    return line.split("!", 1)[0].strip()


def clean_value(value):
    """Clean a scalar namelist value."""
    value = value.strip().rstrip(",")

    if value.startswith('"') and value.endswith('"'):
        return value[1:-1]

    if value.startswith("'") and value.endswith("'"):
        return value[1:-1]

    value_for_float = value.replace("D", "E").replace("d", "e")

    try:
        if any(c in value_for_float for c in [".", "e", "E"]):
            return float(value_for_float)
        return int(value_for_float)
    except ValueError:
        return value


def read_simple_namelist(filename):
    """
    Minimal parser for input.nml.

    It supports lines of the form

        key = value

    and ignores comments after !.
    """
    filename = Path(filename)

    if not filename.exists():
        raise FileNotFoundError(filename)

    cfg = {}
    inside = False

    with filename.open("r") as f:
        for line in f:
            line = strip_comment(line)

            if line == "":
                continue

            if line.startswith("&"):
                inside = True
                continue

            if line.startswith("/"):
                inside = False
                continue

            if not inside:
                continue

            if "=" not in line:
                continue

            key, value = line.split("=", 1)
            key = key.strip().lower()
            value = clean_value(value)

            cfg[key] = value

    return cfg


def get(cfg, key, default):
    return cfg.get(key.lower(), default)


def main():
    input_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("input.nml")

    cfg = read_simple_namelist(input_file)

    nx = int(get(cfg, "nx", 32))
    ny = int(get(cfg, "ny", 32))
    nz = int(get(cfg, "nz", 32))

    xmin = float(get(cfg, "xmin", 0.0))
    xmax = float(get(cfg, "xmax", 1.0))
    ymin = float(get(cfg, "ymin", 0.0))
    ymax = float(get(cfg, "ymax", 1.0))
    zmin = float(get(cfg, "zmin", 0.0))
    zmax = float(get(cfg, "zmax", 1.0))

    bx0 = float(get(cfg, "bx0", 0.0))
    by0 = float(get(cfg, "by0", 0.0))
    bz0 = float(get(cfg, "bz0", 1.0))

    ex0 = float(get(cfg, "ex0", 0.1))
    ey0 = float(get(cfg, "ey0", 0.0))
    ez0 = float(get(cfg, "ez0", 0.0))

    field_file = Path(str(get(cfg, "field_file", "fields_uniform.txt")))

    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    z = np.linspace(zmin, zmax, nz)

    with field_file.open("w") as f:
        f.write("# x y z Bx By Bz Ex Ey Ez\n")
        f.write("# Generated from input.nml by scripts/make_uniform_txt.py\n")
        f.write(f"# nx ny nz = {nx} {ny} {nz}\n")
        f.write(f"# B = ({bx0}, {by0}, {bz0})\n")
        f.write(f"# E = ({ex0}, {ey0}, {ez0})\n")

        # Important: this order matches the Fortran TXT reader:
        # i fastest, then j, then k.
        for zk in z:
            for yj in y:
                for xi in x:
                    f.write(
                        f"{xi:.16e} {yj:.16e} {zk:.16e} "
                        f"{bx0:.16e} {by0:.16e} {bz0:.16e} "
                        f"{ex0:.16e} {ey0:.16e} {ez0:.16e}\n"
                    )

    print(f"Wrote {field_file}")
    print(f"Grid: nx={nx}, ny={ny}, nz={nz}")
    print(f"B = ({bx0}, {by0}, {bz0})")
    print(f"E = ({ex0}, {ey0}, {ez0})")
    print("Columns: x y z Bx By Bz Ex Ey Ez")


if __name__ == "__main__":
    main()

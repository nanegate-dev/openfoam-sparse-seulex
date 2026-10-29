#!/bin/bash
# Minimal reproducible benchmark for OpenFOAM MR !935
# Compares standard dense seulex (sparse false; autonomous false) against
# sparse seulex (sparse true; autonomous true) on OpenFOAM's built-in
# $FOAM_TUTORIALS/combustion/chemFoam tutorials (h2, gri, nc7h16, ic8h18).
#
# Usage:
#   bash reproduce_chemFoam.sh [work_dir] [tutorials...]
# Example:
#   bash reproduce_chemFoam.sh /tmp/seulex_bench h2 gri nc7h16 ic8h18

set -e

WORK_DIR=${1:-/tmp/seulex_mr935_bench}
shift || true
CASES=${*:-h2 gri nc7h16 ic8h18}

mkdir -p "$WORK_DIR"
WORK_DIR=$(cd "$WORK_DIR" && pwd)
echo "Working directory: $WORK_DIR"
echo "Cases: $CASES"
echo "----------------------------------------------------------------------"

for m in $CASES; do
    TUT="$FOAM_TUTORIALS/combustion/chemFoam/$m"
    if [ ! -d "$TUT" ]; then
        echo "Tutorial not found: $TUT (make sure OpenFOAM environment is sourced)"
        exit 1
    fi

    for mode in dense sparse; do
        CASE_DIR="$WORK_DIR/${m}_${mode}"
        rm -rf "$CASE_DIR"
        cp -r "$TUT" "$CASE_DIR"
        cd "$CASE_DIR"

        # Convert CHEMKIN files once if required by the tutorial
        if [ -f chemkin/chem.inp ] && grep -q chemkinToFoam Allrun 2>/dev/null; then
            chemkinToFoam \
                chemkin/chem.inp \
                chemkin/therm.dat \
                chemkin/transportProperties \
                constant/reactions \
                constant/thermo > log.chemkinToFoam 2>&1
        fi

        foamDictionary constant/chemistryProperties -entry odeCoeffs/solver -set seulex > /dev/null
        if [ "$mode" = "dense" ]; then
            foamDictionary constant/chemistryProperties -entry odeCoeffs/sparse -set false > /dev/null
            foamDictionary constant/chemistryProperties -entry odeCoeffs/autonomous -set false > /dev/null
        else
            foamDictionary constant/chemistryProperties -entry odeCoeffs/sparse -set true > /dev/null
            foamDictionary constant/chemistryProperties -entry odeCoeffs/autonomous -set true > /dev/null
        fi

        chemFoam > log.chemFoam 2>&1
        t=$(grep ExecutionTime log.chemFoam | tail -1 | awk '{print $3}')
        printf "%-8s | %-6s : %8s s\n" "$m" "$mode" "$t"
    done

    python3 - "$WORK_DIR/${m}_dense/chemFoam.out" "$WORK_DIR/${m}_sparse/chemFoam.out" <<'PY'
import sys

def read_out(path):
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and not line.startswith("#"):
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
    return rows

a, b = read_out(sys.argv[1]), read_out(sys.argv[2])

def t_ign(rows, Tc=1500.0):
    for (t0, T0), (t1, T1) in zip(rows, rows[1:]):
        if T0 <= Tc < T1:
            return t0 + (Tc - T0) * (t1 - t0) / (T1 - T0)
    return float("nan")

same_steps = (
    len(a) == len(b)
    and all(abs(x[0] - y[0]) <= 1e-12 * max(1.0, abs(x[0])) for x, y in zip(a, b))
)
if same_steps:
    max_dT = max(abs(x[1] - y[1]) for x, y in zip(a, b))
    print(f"   -> Identical {len(a)} steps | max |dT| = {max_dT:.3e} K")
else:
    print(f"   -> Steps: {len(a)} (dense) vs {len(b)} (sparse)")
print(
    f"   -> Ignition time (1500 K): {t_ign(a):.7e} s vs {t_ign(b):.7e} s | "
    f"Final T: {a[-1][1]:.4f} K vs {b[-1][1]:.4f} K"
)
print("----------------------------------------------------------------------")
PY
done

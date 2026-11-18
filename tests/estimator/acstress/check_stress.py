#! /usr/bin/env python3

# Integration check for the ACStress (zero-variance zero-bias stress) estimator.
#
# The same strained diamond system (tests/solids/diamondC_1x1x1_pp_strained_cplx) is run for a few
# VMC blocks with a fixed seed, once with driver_version=legacy and once with driver_version=batch.
# The final block of vmc.s000.scalar.dat is compared with reference values generated from these inputs.
# The two drivers use different random number streams, so each has its own reference table.
# Their agreement with each other was checked statistically, and the multi-walker path is checked
# against the single-walker one in the unit test "ZVZB mw_evaluateStrainDerivsFast matches single walker".
#
# Usage: check_stress.py <legacy|batch> [scalar_file]

import math
import sys

reference_vals = {
    "legacy": {
        "LocalEnergy": -11.633928165,
        "ACStress_hf_0_0": 2.8713954273,
        "ACStress_pulay_0_0": 0.0,
        "ACStress_Ewfgrad_0_0": -5.3698687109,
        "ACStress_wfgrad_0_0": 0.46424324878,
        "ACStress_hf_1_1": -3.8382791532,
        "ACStress_pulay_1_1": 0.0,
        "ACStress_Ewfgrad_1_1": -4.9546419021,
        "ACStress_wfgrad_1_1": 0.41968418624,
        "ACStress_hf_2_2": 1.9324541482,
        "ACStress_pulay_2_2": 0.0,
        "ACStress_Ewfgrad_2_2": -4.3241765938,
        "ACStress_wfgrad_2_2": 0.37282800922,
        "ACStress_hf_0_1": 1.8402576815,
        "ACStress_pulay_0_1": 0.0,
        "ACStress_Ewfgrad_0_1": -0.24023087681,
        "ACStress_wfgrad_0_1": 0.020938579958,
        "ACStress_hf_0_2": -3.5683464602,
        "ACStress_pulay_0_2": 0.0,
        "ACStress_Ewfgrad_0_2": -0.35840878136,
        "ACStress_wfgrad_0_2": 0.026290982325,
        "ACStress_hf_1_2": 0.8055157459,
        "ACStress_pulay_1_2": 0.0,
        "ACStress_Ewfgrad_1_2": 0.33535537366,
        "ACStress_wfgrad_1_2": -0.028699095805,
    },
    "batch": {
        "LocalEnergy": -11.541667959,
        "ACStress_hf_0_0": 2.63845903,
        "ACStress_pulay_0_0": 0.0,
        "ACStress_Ewfgrad_0_0": -4.9071878142,
        "ACStress_wfgrad_0_0": 0.42516643358,
        "ACStress_hf_1_1": -3.7895707904,
        "ACStress_pulay_1_1": 0.0,
        "ACStress_Ewfgrad_1_1": -5.5479415615,
        "ACStress_wfgrad_1_1": 0.4734169906,
        "ACStress_hf_2_2": 1.5227115947,
        "ACStress_pulay_2_2": 0.0,
        "ACStress_Ewfgrad_2_2": -4.8052595614,
        "ACStress_wfgrad_2_2": 0.4165196279,
        "ACStress_hf_0_1": 0.48883984404,
        "ACStress_pulay_0_1": 0.0,
        "ACStress_Ewfgrad_0_1": -0.44698922179,
        "ACStress_wfgrad_0_1": 0.038233908606,
        "ACStress_hf_0_2": -3.7556952789,
        "ACStress_pulay_0_2": 0.0,
        "ACStress_Ewfgrad_0_2": -0.080948965081,
        "ACStress_wfgrad_0_2": 0.0048940940287,
        "ACStress_hf_1_2": 2.455254072,
        "ACStress_pulay_1_2": 0.0,
        "ACStress_Ewfgrad_1_2": 0.28461880271,
        "ACStress_wfgrad_1_2": -0.025946621505,
    },
}

# Relative tolerance for comparing floating point values, plus a tiny absolute one for exact zeros.
REL_TOL = 1e-5
ABS_TOL = 1e-8


def read_last_row(fname):
    with open(fname, "r") as f:
        header = f.readline().lstrip("#").split()
        last = None
        for line in f:
            if line.strip():
                last = line.split()
    if last is None:
        raise RuntimeError("no data rows in " + fname)
    return dict(zip(header, map(float, last)))


if __name__ == "__main__":
    driver = sys.argv[1]
    fname = sys.argv[2] if len(sys.argv) > 2 else "vmc.s000.scalar.dat"
    ref = reference_vals[driver]
    row = read_last_row(fname)

    all_pass = True
    n_stress = 0
    for key, refval in ref.items():
        if key not in row:
            print("Error. Column %s missing from %s" % (key, fname))
            all_pass = False
            continue
        if key.startswith("ACStress"):
            n_stress += 1
        val = row[key]
        if not math.isfinite(val) or not math.isclose(refval, val, rel_tol=REL_TOL, abs_tol=ABS_TOL):
            print("Error. %s Ref = %r Val = %r" % (key, refval, val))
            all_pass = False
        if key.startswith("ACStress_pulay") and val != 0.0:
            print("Error. %s should be exactly zero, got %r" % (key, val))
            all_pass = False

    # 6 independent tensor components x (hf, pulay, Ewfgrad, wfgrad)
    if n_stress != 24:
        print("Error. Expected 24 ACStress columns, found %d" % n_stress)
        all_pass = False

    sys.exit(0 if all_pass else 1)

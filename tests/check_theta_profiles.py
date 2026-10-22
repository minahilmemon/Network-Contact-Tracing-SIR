#!/usr/bin/env python3
"""Run: python3 tests/check_theta_profiles.py

Checks the analytic no-tracing solution, direct quadrature of the original
pair characteristics, probability bounds, and age/degree/tolerance refinement.
Writes reproducible numerical comparisons beside the default profile outputs.
"""
from dataclasses import replace
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from theta_profiles import Parameters, diagnostics, no_tracing_theta, solve_profiles


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    p = Parameters()
    reference = solve_profiles(replace(p, p=0), age_max=100, tolerance=1e-9)
    ages = np.linspace(0, 100, 1001)
    analytic_error = float(np.max(abs(reference.values(ages)[:p.max_degree]
                                     - no_tracing_theta(p, ages))))
    rate_error = abs(reference.growth_rate-(p.beta*p.kappa-p.beta-p.alpha-p.gamma))
    check(analytic_error < 1e-8 and rate_error < 1e-8,
          "No-tracing profiles or growth rate disagree with the analytic reference.")
    print(f"PASS analytic p=0: profile error {analytic_error:.3g}, "
          f"growth-rate error {rate_error:.3g}", flush=True)

    # With alpha=0 and no diagnosed seed, tracing has no source, even if p>0.
    no_diagnosis = solve_profiles(replace(p, alpha=0), tolerance=1e-8)
    check(np.max(abs(no_diagnosis.values(ages)[:p.max_degree])) < 1e-9,
          "Tracing appeared without a diagnosis source.")
    check(abs(no_diagnosis.growth_rate-(p.beta*p.kappa-p.beta-p.gamma)) < 1e-8,
          "Incorrect growth rate without diagnosis.")
    print("PASS no diagnosis: identically zero tracing profile.", flush=True)

    base = solve_profiles(p, tolerance=1e-8)
    report = diagnostics(base)
    for key in ["normalization_error", "scalar_identity_error", "max_initial_slope_error",
                "max_direct_theta_equation_error", "max_pair_transport_error"]:
        check(report[key] < 1e-8, f"Original-equation diagnostic failed: {key}")
    check(report["max_relative_edge_symmetry_error"] < 1e-10,
          "Reconstructed edge densities are not symmetric.")
    check(report["min_theta"] >= -1e-10 and report["max_theta"] <= 1,
          "Tracing probabilities outside [0,1].")
    check(report["max_grid_S_I_T_probability_mass"] < 1+1e-7,
          "Neighbour probability mass exceeds one.")
    check(report["endpoint_R_flux"] < 1e-12, "Infected-neighbour tail is too large.")
    print(f"PASS original equations, normalization, probability bounds, and edge symmetry; "
          f"max Theta residual {report['max_direct_theta_equation_error']:.3g}", flush=True)

    comparisons = []
    configurations = [
        ("age cutoff 60", p, dict(age_max=60), 2e-6),
        ("age cutoff 140", p, dict(age_max=140), 1e-8),
        ("tolerance 1e-9", p, dict(tolerance=1e-9), 1e-8),
        ("initial spacing 0.2", p, dict(initial_step=.2), 1e-8),
        ("single tracing continuation step", p, dict(continuation_step=1), 1e-8),
        ("degree cutoff 30", replace(p, max_degree=30), {}, 2e-6),
    ]
    for label, parameters, overrides, allowed in configurations:
        options = dict(tolerance=1e-8)
        options.update(overrides)
        other = solve_profiles(parameters, **options)
        shared_ages = np.linspace(0, min(base.age_max, other.age_max), 1001)
        profile_change = float(np.max(abs(
            base.values(shared_ages)[:p.max_degree]
            - other.values(shared_ages)[:p.max_degree])))
        rate_change = abs(base.growth_rate-other.growth_rate)
        check(profile_change < allowed and rate_change < allowed,
              f"Refinement failed: {label}, profile={profile_change:g}, rate={rate_change:g}")
        comparisons.append(dict(case=label, growth_rate=other.growth_rate,
                                max_theta_change=profile_change, growth_rate_change=rate_change,
                                nodes=len(other.result.x)))
        print(f"PASS {label}: max profile change {profile_change:.3g}, "
              f"growth-rate change {rate_change:.3g}", flush=True)
    output = ROOT/"output"/"theta-profiles"/"verification.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(dict(analytic_profile_error=analytic_error,
                                     analytic_growth_rate_error=rate_error,
                                     original_equation_diagnostics=report,
                                     refinement=comparisons), indent=2)+"\n")
    print(f"All checks passed. {output}")


if __name__ == "__main__":
    main()

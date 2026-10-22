#!/usr/bin/env python3
"""Compare the stationary-profile growth rate with the full nonlinear PDE.

First run theta_profiles.py, then tests/check_theta_full_pde.R at
--dt=0.4, 0.2, 0.1, and 0.05. This script reads those outputs and the
existing finite-seed epidemic CSV; it never changes the epidemic inputs.
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/"output"/"theta-profiles"/"full-pde"


def read_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def main():
    reference = json.loads((OUTPUT.parent/"diagnostics.json").read_text())
    target = reference["growth_rate"]
    rows = []
    for step in [.4, .2, .1, .05]:
        fits = read_csv(OUTPUT/f"growth-dt{step:g}.csv")
        fit = fits[(fits["start"] == 25) & (fits["end"] == 30)][0]
        rows.append(dict(dt=step, growth_Lambda=float(fit["growth_Lambda"]),
                         growth_I=float(fit["growth_I"]),
                         absolute_error=abs(float(fit["growth_Lambda"])-target),
                         susceptible_loss=float(fit["susceptible_loss"]),
                         max_log_I_fit_residual=float(fit["max_log_I_fit_residual"])))
    extrapolations = []
    for coarse, fine in zip(rows[:-1], rows[1:]):
        rate = 2*fine["growth_Lambda"]-coarse["growth_Lambda"]
        extrapolations.append(dict(coarse_dt=coarse["dt"], fine_dt=fine["dt"],
                                   extrapolated_growth_rate=rate, absolute_error=abs(rate-target)))
    # The PDE method is first order; the checks require decreasing errors
    # and agreement after eliminating its leading time/age-grid error.
    errors = np.array([row["absolute_error"] for row in rows])
    if not np.all(np.diff(errors) < 0):
        raise AssertionError("PDE growth does not converge toward the stationary-profile prediction.")
    if extrapolations[-1]["absolute_error"] > 2e-5:
        raise AssertionError("Refined extrapolation disagrees with the predicted growth rate.")
    if max(row["susceptible_loss"] for row in rows) > 1e-7:
        raise AssertionError("Susceptible depletion is too large for an invasion-growth comparison.")
    if max(abs(row["growth_I"]-row["growth_Lambda"]) for row in rows) > 3e-6:
        raise AssertionError("Node and exposure growth rates have not settled to the same mode.")

    saved = read_csv(ROOT/"PDE_contact_tracing_dt0.1.csv")
    finite_seed_fits = []
    for lo, hi in [(0,5), (5,10), (10,15), (15,20), (20,25)]:
        mask = (saved["t"] >= lo) & (saved["t"] <= hi)
        rate = np.polyfit(saved["t"][mask], np.log(saved["I_PDE"][mask]), 1)[0]
        finite_seed_fits.append(dict(start=lo, end=hi, growth_I=float(rate),
                                    final_susceptible_fraction=float(saved["S_PDE"][mask][-1]/100000)))
    report = dict(predicted_growth_rate=target, parameters=reference["parameters"],
                  initial_infected_fraction=1e-11, fitting_window=[25,30],
                  fitted_quantity="log of mean Lambda across degree classes",
                  full_pde=rows, first_order_extrapolations=extrapolations,
                  saved_finite_seed_trajectory=finite_seed_fits)
    (OUTPUT/"growth-validation.json").write_text(json.dumps(report, indent=2)+"\n")
    with (OUTPUT/"growth-refinement.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    ax = axes[0]
    tiny = read_csv(OUTPUT/"trajectory-dt0.1.csv")
    for trajectory, label, colour in [
        (saved, "Saved epidemic: I(0)/N = 0.0002", "#bf6b18"),
        (tiny, "Tiny seed: I(0)/N = 10⁻¹¹", "#2679a8"),
    ]:
        slope = np.gradient(np.log(trajectory["I_PDE"]), trajectory["t"])
        ax.plot(trajectory["t"], slope, label=label, color=colour, linewidth=2)
    ax.axhline(target, color="#333333", linestyle="--", label="Stationary-profile prediction")
    ax.set(xlim=(5,30), ylim=(.095,.125), xlabel="Time",
           ylabel="Instantaneous growth rate of I", title="Full PDE with dt = da = 0.1")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.grid(alpha=.2)

    ax = axes[1]
    steps = np.array([row["dt"] for row in rows])
    rates = np.array([row["growth_Lambda"] for row in rows])
    ax.plot(steps, rates, "o-", color="#2679a8", label="Full PDE: tiny seed")
    ax.axhline(target, color="#333333", linestyle="--",
               label=f"Profile prediction: {target:.8f}")
    extrapolated = extrapolations[-1]["extrapolated_growth_rate"]
    ax.plot([0,steps[-1]], [extrapolated,rates[-1]], ":", color="#2679a8")
    ax.plot(0, extrapolated, "D", color="#2679a8", markersize=5,
            label=f"First-order extrapolation: {extrapolated:.8f}")
    ax.set(xlim=(-.015,.415), xlabel="Time and age step (dt = da)",
           ylabel="Fitted growth rate, times 25–30", title="Convergence toward the predicted rate")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.grid(alpha=.2)
    fig.savefig(OUTPUT/"growth-comparison.png", dpi=180)
    plt.close(fig)
    print(json.dumps(report, indent=2))
    print("PASS full-PDE growth converges to the stationary-profile prediction.")


if __name__ == "__main__":
    main()

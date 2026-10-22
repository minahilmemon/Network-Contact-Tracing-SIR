#!/usr/bin/env python3
"""Stationary invasion-limit conditional profiles for uncorrelated mixing.

Run: python3 theta_profiles.py
Dependencies: numpy, scipy, matplotlib. No epidemic simulation is run.

The unknowns are Theta_k(a), T_k(a)=integral_0^a Theta_k, two auxiliary
integrals U,V, an accumulated normalization J, and a scalar s. Writing
g=delta*p, h_j=g*(j-1), mu=alpha+gamma, and w_j=alpha+h_j*Theta_j,

 U(a) = sum_j q_j integral_0^a w_j(b) exp[-beta(a-b)-mu*b-h_j*T_j(b)] db
 V(a) = sum_j (j-1)*q_j exp[-mu*a]
                     integral_a^L w_j(b) exp[s*(b-a)-h_j*T_j(b)] db.

Differentiating these expressions gives a local boundary-value system.
Theta_k(0)=T_k(0)=U(0)=J(0)=0, V(L)=0, J(L)=1. The normalizing scalar
s is solved with the profiles; lambda=-s-beta-alpha-gamma is computed
afterward. R_j=beta*(j-1)*q_j*exp[s*a-h_j*T_j]. Degree zero is excluded
because these probabilities condition on a contact.

Age truncation L replaces infinite-age integrals by integrals to L. In
particular, s+beta*kappa=A+sum_j R_j(L), where A=sum_j integral h_j*Theta_j*R_j.
The endpoint flux and cutoff refinement measure this approximation.
Continuation from p=0 selects one numerical branch; convergence is not
an existence or uniqueness proof. See tests/check_theta_profiles.py.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, cumulative_trapezoid, quad_vec, solve_bvp
from scipy.stats import poisson


@dataclass(frozen=True)
class Parameters:
    beta: float = 0.10
    alpha: float = 0.15
    gamma: float = 0.10
    delta: float = 1.0
    p: float = 0.3
    mean_degree: float = 5.0
    max_degree: int = 20

    def __post_init__(self):
        rates = np.array([self.beta, self.alpha, self.gamma, self.delta,
                          self.p, self.mean_degree])
        if (not np.all(np.isfinite(rates)) or self.beta <= 0
                or min(self.alpha, self.gamma, self.delta) < 0
                or not 0 <= self.p <= 1 or self.mean_degree <= 0
                or self.max_degree < 2):
            raise ValueError("Require beta, mean_degree > 0, nonnegative rates, "
                             "0 <= p <= 1, and max_degree >= 2.")

    @property
    def degrees(self):
        return np.arange(1, self.max_degree + 1)

    @property
    def q(self):
        weights = self.degrees * poisson.pmf(self.degrees, self.mean_degree)
        return weights / weights.sum()

    @property
    def kappa(self):
        return float(self.q @ (self.degrees - 1))


def decay_convolution(first, second, age):
    """Integral_0^a exp[-first*b-second*(a-b)] db, including equal rates."""
    difference = abs(first - second)
    if difference < 1e-12:
        return age * np.exp(-first * age)
    return np.exp(-min(first, second) * age) * (-np.expm1(-difference * age)) / difference


def no_tracing_theta(parameters, age):
    """Infinite-age analytic p=0 reference (diagnosis remains active)."""
    p = parameters
    mu = p.alpha + p.gamma
    if abs(mu - p.beta) > 1e-8:
        double_convolution = (decay_convolution(p.beta, p.delta, age)
                              - decay_convolution(mu, p.delta, age)) / (mu - p.beta)
    else:
        # A repeated exponential rate; quadrature avoids cancellation.
        double_convolution = np.array([
            quad_vec(lambda b: b * np.exp(-p.beta*b-p.delta*(a-b)), 0, a)[0]
            for a in age])
    k = p.degrees[:, None]
    return (p.alpha / k * decay_convolution(mu, p.delta, age)
            + p.beta * (k-1) / k * p.alpha * double_convolution)


def equations(parameters):
    p = parameters
    K = p.max_degree
    k = p.degrees[:, None]
    q = p.q
    excess_q = q * (p.degrees-1)
    g = p.delta*p.p
    h = g*(k-1)
    mu = p.alpha+p.gamma

    def rhs(age, values, scalar):
        theta, integral = values[:K], values[K:2*K]
        U, V = values[-3], values[-2]
        neighbour_source = (p.alpha+h*theta) * np.exp(-mu*age-h*integral)
        theta_prime = (-p.delta*theta+g*theta**2
                       + p.beta/k*np.exp(g*integral)*((k-1)*U+V))
        return np.vstack((theta_prime, theta,
                          -p.beta*U+q@neighbour_source,
                          -(mu+scalar[0])*V-excess_q@neighbour_source,
                          p.beta*excess_q@np.exp(scalar[0]*age-h*integral)))

    def boundary(left, right, scalar):
        return np.r_[left[:2*K], left[-3], left[-1], right[-2], right[-1]-1]

    return rhs, boundary


@dataclass
class ProfileSolution:
    parameters: Parameters
    age_max: float
    tolerance: float
    result: object
    continuation: list

    @property
    def s(self):
        return float(self.result.p[0])

    @property
    def growth_rate(self):
        p = self.parameters
        return -self.s-p.beta-p.alpha-p.gamma

    def values(self, age):
        ages = np.asarray(age)
        if np.any(ages < 0) or np.any(ages > self.age_max):
            raise ValueError("Age outside the solved interval.")
        return self.result.sol(age)

    def R(self, age):
        p = self.parameters
        age = np.atleast_1d(age)
        integral = self.values(age)[p.max_degree:2*p.max_degree]
        h = p.delta*p.p*(p.degrees-1)
        return (p.beta*(p.degrees-1)*p.q)[:, None] * np.exp(
            self.s*age-h[:, None]*integral)

    def pairs(self, k, focal_age, neighbour_age):
        """All j at fixed k,a,b, reconstructed from the original characteristics.

        No pair array is retained by the solver. At the corner/diagonal the
        focal-infection boundary is used; its point value does not affect
        any age integral. The two corner inflows need not agree.
        """
        p = self.parameters
        if not 1 <= k <= p.max_degree:
            raise ValueError("Degree outside the solved classes.")
        a, b = focal_age, neighbour_age
        g = p.delta*p.p
        h = g*(p.degrees-1)
        K = p.max_degree
        Ta = self.values(a)[K+k-1]
        Tb = self.values(b)[K:2*K]
        mu = p.alpha+p.gamma
        if b >= a:
            shifted = self.values(b-a)[K:2*K]
            return self.R(b-a)[:, 0]/k * np.exp(-mu*a-h*(Tb-shifted)+g*Ta)
        shifted = self.values(a-b)[K+k-1]
        C = (k-1)/k*p.q*np.exp(-p.beta*(a-b)+g*shifted)
        return p.beta*C*np.exp(-mu*b-h*Tb+g*(Ta-shifted))


def solve_profiles(parameters=Parameters(), age_max=100.0, tolerance=1e-8,
                   initial_step=0.5, max_nodes=20000, continuation_step=0.1,
                   progress=False):
    if (not np.isfinite(age_max) or age_max <= 0
            or not np.isfinite(tolerance) or tolerance <= 0
            or not np.isfinite(initial_step) or initial_step <= 0
            or not np.isfinite(continuation_step) or continuation_step <= 0):
        raise ValueError("Age cutoff, tolerance, and step sizes must be positive.")
    p = parameters
    age = np.linspace(0, age_max, max(3, int(np.ceil(age_max/initial_step))+1))
    theta = no_tracing_theta(p, age)
    integral = cumulative_trapezoid(theta, age, axis=1, initial=0)
    mu = p.alpha+p.gamma
    U = p.alpha*decay_convolution(p.beta, mu, age)
    V = p.alpha/p.beta*np.exp(-mu*age)
    scalar = np.array([-p.beta*p.kappa])
    values = np.vstack((theta, integral, U, V, -np.expm1(scalar[0]*age)))
    trace_probabilities = np.linspace(0, p.p, max(1, int(np.ceil(p.p/continuation_step)))+1)
    trace_probabilities = np.unique(trace_probabilities)
    history = []
    for probability in trace_probabilities:
        stage = replace(p, p=float(probability))
        rhs, boundary = equations(stage)
        result = solve_bvp(rhs, boundary, age, values, p=scalar,
                           tol=tolerance, bc_tol=tolerance, max_nodes=max_nodes)
        if not result.success:
            raise RuntimeError(f"BVP failed at p={probability:g}: {result.message}")
        if not np.all(np.isfinite(result.y)) or not np.all(np.isfinite(result.p)):
            raise RuntimeError("BVP returned nonfinite values.")
        history.append(dict(p=float(probability), nodes=len(result.x),
                            max_scaled_residual=float(max(result.rms_residuals)),
                            s=float(result.p[0])))
        if progress:
            print(f"p={probability:.3f}: {len(result.x)} nodes, "
                  f"residual={max(result.rms_residuals):.2g}, s={result.p[0]:.10f}", flush=True)
        age, values, scalar = result.x, result.y, result.p
    return ProfileSolution(p, age_max, tolerance, result, history)


def diagnostics(solution):
    """Independent quadrature of the characteristic solution, not U/V ODEs."""
    p = solution.parameters
    K = p.max_degree
    g = p.delta*p.p
    h = g*(p.degrees-1)
    L = solution.age_max
    age = np.linspace(0, L, max(2001, int(np.ceil(L/.025))+1))
    values = solution.values(age)
    theta = values[:K]
    integral = values[K:2*K]
    # Whole-grid probability mass from separate quadrature with w_j=1.
    # Reverse accumulation avoids cancellation in tiny old-age tails.
    survival = np.exp(-(p.alpha+p.gamma)*age-h[:, None]*integral)
    U_mass = np.exp(-p.beta*age)*cumulative_simpson(
        np.exp(p.beta*age)*(p.q@survival), x=age, initial=0)
    weighted_R = (p.q*(p.degrees-1)) @ np.exp(solution.s*age-h[:, None]*integral)
    V_mass = np.exp(-(p.alpha+p.gamma+solution.s)*age)*cumulative_simpson(
        weighted_R[::-1], x=-age[::-1], initial=0)[::-1]
    k_grid = p.degrees[:, None]
    B_mass = p.beta/k_grid*np.exp(g*integral)*((k_grid-1)*U_mass+V_mass)
    C_mass = (k_grid-1)/k_grid*np.exp(-p.beta*age+g*integral)
    total_mass = theta+B_mass+C_mass
    mass = float(quad_vec(lambda b: solution.R(b)[:, 0].sum(), 0, L,
                          epsabs=1e-10, epsrel=1e-10)[0])
    A = float(quad_vec(lambda b: (h*solution.values(b)[:K]*solution.R(b)[:, 0]).sum(),
                       0, L, epsabs=1e-10, epsrel=1e-10)[0])
    endpoint_flux = float(solution.R(L).sum())
    source_errors, probability_masses, slopes, pair_errors, edge_errors = [], [], [], [], []
    degrees = sorted(set([1, min(2, K), min(5, K), K]))
    for k in degrees:
        for a in [0., 1., 3., 10., min(30., L/2), min(60., .8*L)]:
            if a > L:
                continue

            def integrand(b):
                B = solution.pairs(k, a, b)
                return np.array([((p.alpha+h*solution.values(b)[:K])*B).sum(), B.sum()])

            integrated = quad_vec(integrand, 0, a, epsabs=1e-10, epsrel=1e-9)[0]
            integrated += quad_vec(integrand, a, L, epsabs=1e-10, epsrel=1e-9)[0]
            value = solution.values(a)
            th = value[k-1]
            rhs = -p.delta*th+g*th**2+integrated[0]
            source_errors.append(abs(solution.result.sol(a, 1)[k-1]-rhs))
            C = (k-1)/k*np.exp(-p.beta*a+g*value[K+k-1])
            probability_masses.append(C+th+integrated[1])
        slopes.append(abs(solution.result.sol(0, 1)[k-1]-(p.alpha*mass+A)/k))
        for a, b in [(1., 3.), (4., 1.)]:
            if max(a, b)+.001 > L:
                continue
            eps = 1e-4
            derivative = ((solution.pairs(k, a+eps, b+eps)
                           - solution.pairs(k, a-eps, b-eps))/(2*eps))
            B = solution.pairs(k, a, b)
            rhs = (-(p.alpha+p.gamma)-h*solution.values(b)[:K]
                   +g*solution.values(a)[k-1])*B
            pair_errors.append(float(np.max(abs(derivative-rhs))))
            # Symmetry of the underlying edge density: k*i_k(a)*B_kj(a,b)
            # equals j*i_j(b)*B_jk(b,a), with i_k proportional to k*p_k.
            for j in degrees:
                node_k = k*poisson.pmf(k, p.mean_degree)*np.exp(
                    -(solution.growth_rate+p.alpha+p.gamma)*a
                    -g*k*solution.values(a)[K+k-1])
                node_j = j*poisson.pmf(j, p.mean_degree)*np.exp(
                    -(solution.growth_rate+p.alpha+p.gamma)*b
                    -g*j*solution.values(b)[K+j-1])
                left = k*node_k*B[j-1]
                right = j*node_j*solution.pairs(j, b, a)[k-1]
                edge_errors.append(abs(left-right)/max(abs(left), abs(right), 1e-30))
    peaks = []
    for k in p.degrees:
        index = int(np.argmax(theta[k-1]))
        peaks.append(dict(degree=int(k), peak=float(theta[k-1, index]),
                          peak_age=float(age[index])))
    return dict(parameters=asdict(p), age_max=L, tolerance=solution.tolerance,
                s=solution.s, growth_rate=solution.growth_rate, mean_tracing_hazard=A,
                kappa=p.kappa, omitted_neighbour_degree_mass=float(poisson.sf(K-1, p.mean_degree)),
                normalization=mass, normalization_error=abs(mass-1),
                endpoint_R_flux=endpoint_flux,
                scalar_identity_error=abs(solution.s+p.beta*p.kappa-A-endpoint_flux),
                max_scaled_collocation_residual=float(max(solution.result.rms_residuals)),
                max_direct_theta_equation_error=float(max(source_errors)),
                max_initial_slope_error=float(max(slopes)),
                max_pair_transport_error=float(max(pair_errors, default=0)),
                max_relative_edge_symmetry_error=float(max(edge_errors, default=0)),
                min_theta=float(theta.min()), max_theta=float(theta.max()),
                max_grid_S_I_T_probability_mass=float(total_mass.max()),
                max_sampled_S_I_T_probability_mass=float(max(probability_masses)),
                theta_peaks=peaks, continuation=solution.continuation)


def save_outputs(solution, report, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    age = np.linspace(0, solution.age_max, int(np.ceil(solution.age_max/.1))+1)
    K = solution.parameters.max_degree
    values = solution.values(age)
    R = solution.R(age)
    with (output_dir/"profiles.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["infection_age", "degree", "theta", "integrated_theta", "R"])
        for index, a in enumerate(age):
            for k in range(1, K+1):
                writer.writerow([a, k, values[k-1, index], values[K+k-1, index], R[k-1, index]])
    (output_dir/"diagnostics.json").write_text(json.dumps(report, indent=2)+"\n")
    fig, ax = plt.subplots(figsize=(8.2, 4.9), layout="constrained")
    for k in sorted(set([1, min(2, K), min(5, K), min(10, K), K])):
        ax.plot(age, values[k-1], label=f"k = {k}", linewidth=2)
    ax.set(xlabel="Infection age", ylabel=r"Traced-neighbour probability $\Theta_k(a)$",
           xlim=(0, min(60, solution.age_max)), ylim=(0, None),
           title="Stationary conditional profiles near the disease-free state")
    ax.legend(title="Focal degree", frameon=False)
    ax.grid(alpha=.18)
    p = solution.parameters
    ax.text(.98, .97, f"Poisson mean degree {p.mean_degree:g}\n"
            f"β={p.beta:g}, α={p.alpha:g}, γ={p.gamma:g}\nδ={p.delta:g}, p={p.p:g}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9)
    fig.savefig(output_dir/"theta-profiles.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ["beta", "alpha", "gamma", "delta", "p", "mean_degree"]:
        parser.add_argument("--"+name.replace("_", "-"), type=float, default=getattr(Parameters(), name))
    parser.add_argument("--max-degree", type=int, default=20)
    parser.add_argument("--age-max", type=float, default=100)
    parser.add_argument("--tol", type=float, default=1e-8)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parent/"output"/"theta-profiles")
    args = parser.parse_args()
    p = Parameters(**{key: getattr(args, key) for key in asdict(Parameters())})
    solution = solve_profiles(p, args.age_max, args.tol, progress=True)
    report = diagnostics(solution)
    # These checks reject obvious nonphysical branches; sampled bounds do not
    # prove admissibility on the whole infinite interval.
    if (report["min_theta"] < -1e-7 or report["max_theta"] > 1+1e-7
            or report["max_sampled_S_I_T_probability_mass"] > 1+1e-6
            or report["max_grid_S_I_T_probability_mass"] > 1+1e-6
            or report["normalization_error"] > 1e-6
            or report["max_direct_theta_equation_error"] > 1e-6):
        raise RuntimeError("Profile diagnostics failed; no outputs written.\n"+json.dumps(report, indent=2))
    save_outputs(solution, report, args.output_dir)
    print(f"Normalization = {report['normalization']:.12f}; "
          f"direct equation error = {report['max_direct_theta_equation_error']:.3g}")
    print(f"Derived growth rate = {solution.growth_rate:.10f}; "
          f"mean tracing hazard = {report['mean_tracing_hazard']:.10f}")
    print(f"Outputs: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()

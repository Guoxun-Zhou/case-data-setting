"""Static switch-contribution baseline adapted from Gao et al. (2024).

The program deliberately reproduces only the *contribution-screening* idea in
Gao et al., IEEE Transactions on Power Systems, 39(6), 7064-7076, rather than
their WQMIX training procedure.  It is therefore a transparent model-based
baseline for comparison with a model-based spatiotemporal reduction method.

Differences from Gao et al. are explicit: their PV-curtailment / load-shedding
indices are replaced by the two objectives used in this study, i.e. network
loss and branch-current load imbalance.  The output is one static set of K
key switches.  This baseline does not divide the day into periods.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import skfuzzy as fuzz
from gurobipy import GRB, Model, quicksum


LOSS_WEIGHT = 1.0
BALANCE_WEIGHT = 7.0
SCALE = 1000.0
BASE_STATUS = np.array([1] * 13 + [0] * 3, dtype=int)


@dataclass
class CaseData:
    branch_in: np.ndarray
    branch_out: np.ndarray
    r: np.ndarray
    x: np.ndarray
    p: np.ndarray
    q: np.ndarray
    pv: np.ndarray
    sources: np.ndarray

    @property
    def n_branch(self) -> int:
        return len(self.r)

    @property
    def n_node(self) -> int:
        return self.p.shape[0]

    @property
    def n_time(self) -> int:
        return self.p.shape[1]


def read_case(case_dir: Path) -> CaseData:
    """Read exactly the 16-node files used by the manuscript."""
    s_base, v_base = 5.68, 12.66
    r_base = v_base**2 / s_base
    branch = pd.read_excel(case_dir / "branch.xlsx", index_col=0).values
    p = pd.read_excel(case_dir / "load_P.xlsx").values / s_base / 1000.0
    q = pd.read_excel(case_dir / "load_Q.xlsx").values[: p.shape[0], :] / s_base / 1000.0
    pv = pd.read_excel(case_dir / "PV.xlsx").values / s_base / 1000.0
    return CaseData(
        branch_in=branch[:, 0].astype(int),
        branch_out=branch[:, 1].astype(int),
        r=branch[:, 2] / r_base,
        x=branch[:, 3] / r_base,
        p=p,
        q=q,
        pv=pv,
        sources=np.arange(3, dtype=int),
    )


def make_samples(data: CaseData, sample_count: int, seed: int) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Return one nominal trace and reproducible modest daily operating variations."""
    rng = np.random.default_rng(seed)
    samples = [(data.p.copy(), data.q.copy(), data.pv.copy())]
    for _ in range(sample_count - 1):
        # Common daily conditions and small node/hour deviations retain the
        # original 24-hour temporal pattern while making M scenarios.
        load_common = rng.uniform(0.90, 1.10)
        pv_common = rng.uniform(0.85, 1.15)
        load_local = rng.uniform(0.95, 1.05, size=data.p.shape)
        pv_local = rng.uniform(0.95, 1.05, size=data.pv.shape)
        samples.append((data.p * load_common * load_local,
                        data.q * load_common * load_local,
                        data.pv * pv_common * pv_local))
    return samples


def solve_daily_opf(data: CaseData, p: np.ndarray, q: np.ndarray, pv: np.ndarray,
                    fixed_status: dict[int, int] | None = None) -> dict[str, np.ndarray | float]:
    """Solve the manuscript's 24-hour OPF with independent hourly topologies.

    ``fixed_status`` fixes a branch for every hour.  Passing BASE_STATUS yields
    the no-reconfiguration reference; passing None yields the full-switch
    reference used to collect Gao-style action samples.
    """
    fixed_status = fixed_status or {}
    b_count, n_count, t_count = data.n_branch, data.n_node, data.n_time
    m = Model("gao_static_contribution")
    m.Params.OutputFlag = 0
    m.Params.NonConvex = 2
    z = m.addVars(b_count, t_count, vtype=GRB.BINARY, name="z")
    fz = m.addVars(b_count, t_count, lb=0.0, ub=1.0, name="fz")
    ff = m.addVars(b_count, t_count, lb=0.0, ub=1.0, name="ff")
    node_p = m.addVars(n_count, t_count, lb=-GRB.INFINITY, name="P")
    node_q = m.addVars(n_count, t_count, lb=-GRB.INFINITY, name="Q")
    load_shed_p = m.addVars(n_count, t_count, lb=-GRB.INFINITY, name="Plr")
    load_shed_q = m.addVars(n_count, t_count, lb=-GRB.INFINITY, name="Qlr")
    pv_output = m.addVars(n_count, t_count, lb=-GRB.INFINITY, name="Ppv")
    voltage_sq = m.addVars(n_count, t_count, lb=0.0, name="U")
    branch_p = m.addVars(b_count, t_count, lb=-GRB.INFINITY, name="BP")
    branch_q = m.addVars(b_count, t_count, lb=-GRB.INFINITY, name="BQ")
    current_sq = m.addVars(b_count, t_count, lb=0.0, name="I")

    for branch, state in fixed_status.items():
        m.addConstrs(z[branch, t] == state for t in range(t_count))

    source_set = set(data.sources.tolist())
    nonsource_count = n_count - len(source_set)
    for t in range(t_count):
        for node in range(n_count):
            m.addConstr(pv_output[node, t] <= 0.0)
            m.addConstr(pv_output[node, t] >= pv[node, t])
            if node in source_set:
                m.addConstr(load_shed_p[node, t] == 0.0)
                m.addConstr(load_shed_q[node, t] == 0.0)
            else:
                m.addConstr(load_shed_p[node, t] >= 0.0)
                m.addConstr(load_shed_p[node, t] <= 0.1 * p[node, t])
                m.addConstr(load_shed_q[node, t] == load_shed_p[node, t] * q[node, t] / p[node, t])
                m.addConstr(node_p[node, t] == p[node, t] - load_shed_p[node, t] + pv_output[node, t])
                m.addConstr(node_q[node, t] == q[node, t] - load_shed_q[node, t])

        m.addConstrs(current_sq[b, t] <= 10.0 * z[b, t] for b in range(b_count))
        m.addConstrs(fz[b, t] + ff[b, t] == z[b, t] for b in range(b_count))
        source_fz = 0.0
        source_ff = 0.0
        for node in range(n_count):
            outgoing = [b for b in range(b_count) if data.branch_in[b] == node]
            incoming = [b for b in range(b_count) if data.branch_out[b] == node]
            if node not in source_set:
                m.addConstr(quicksum(fz[b, t] for b in incoming) + quicksum(ff[b, t] for b in outgoing)
                            <= nonsource_count / (nonsource_count + 1))
            else:
                source_fz += quicksum(fz[b, t] for b in incoming)
                source_ff += quicksum(ff[b, t] for b in outgoing)
            m.addConstr(quicksum(branch_p[b, t] for b in outgoing) + node_p[node, t]
                        == quicksum(branch_p[b, t] - data.r[b] * current_sq[b, t] for b in incoming))
            m.addConstr(quicksum(branch_q[b, t] for b in outgoing) + node_q[node, t]
                        == quicksum(branch_q[b, t] - data.x[b] * current_sq[b, t] for b in incoming))
        m.addConstr(source_fz + source_ff <= nonsource_count / (nonsource_count + 1))

        for b in range(b_count):
            i, j = data.branch_in[b], data.branch_out[b]
            drop = -voltage_sq[j, t] + voltage_sq[i, t] - 2 * (data.r[b] * branch_p[b, t] + data.x[b] * branch_q[b, t]) + (data.r[b]**2 + data.x[b]**2) * current_sq[b, t]
            m.addConstr(drop <= 100.0 * (1 - z[b, t]))
            m.addConstr(drop >= -100.0 * (1 - z[b, t]))
            m.addConstr(4 * branch_p[b, t]**2 + 4 * branch_q[b, t]**2 + (current_sq[b, t] - voltage_sq[i, t])**2
                        <= (current_sq[b, t] + voltage_sq[i, t])**2)
        m.addConstrs(voltage_sq[s, t] == 1.05**2 for s in data.sources)
        m.addConstrs(voltage_sq[n, t] <= 1.05**2 for n in range(n_count))
        m.addConstrs(voltage_sq[n, t] >= 0.95**2 for n in range(n_count))

        for node in data.sources:
            m.addConstr(p[node, t] - 200.0 - node_p[node, t] <= 0.0)
            m.addConstr(q[node, t] - 200.0 - node_q[node, t] <= 0.0)
            m.addConstr(p[node, t] - node_p[node, t] >= 0.0)
            m.addConstr(q[node, t] - node_q[node, t] + 10000.0 >= 0.0)

    loss = quicksum(SCALE * current_sq[b, t] * data.r[b] for b in range(b_count) for t in range(t_count))
    imbalance = quicksum(SCALE * quicksum((current_sq[b, t] - quicksum(current_sq[k, t] for k in range(b_count)) / b_count) ** 2
                                           for b in range(b_count)) / b_count for t in range(t_count))
    pv_penalty = quicksum(pv_output[n, t] - pv[n, t] for n in range(n_count) for t in range(t_count))
    m.setObjective(LOSS_WEIGHT * loss + BALANCE_WEIGHT * imbalance + pv_penalty, GRB.MINIMIZE)
    m.optimize()
    if m.Status not in (GRB.OPTIMAL, GRB.SUBOPTIMAL):
        raise RuntimeError(f"Gurobi did not solve the model; status={m.Status}")
    z_value = np.array([[round(z[b, t].X) for t in range(t_count)] for b in range(b_count)], dtype=int)
    is_value = np.array([[current_sq[b, t].X for t in range(t_count)] for b in range(b_count)])
    loss_t = SCALE * np.sum(is_value * data.r[:, None], axis=0)
    imbalance_t = SCALE * np.mean((is_value - np.mean(is_value, axis=0, keepdims=True)) ** 2, axis=0)
    return {"z": z_value, "loss_t": loss_t, "imbalance_t": imbalance_t,
            "objective_t": LOSS_WEIGHT * loss_t + BALANCE_WEIGHT * imbalance_t,
            "objective": float(m.ObjVal)}


def contribution_scores(data: CaseData, sample_count: int, seed: int) -> tuple[np.ndarray, list[dict[str, object]]]:
    """Calculate the adapted Gao-style static score (Eq. 19--23 analogue)."""
    scores = np.zeros(data.n_branch)
    trace: list[dict[str, object]] = []
    base_fixed = {b: int(BASE_STATUS[b]) for b in range(data.n_branch)}
    for sample_id, (p, q, pv) in enumerate(make_samples(data, sample_count, seed)):
        baseline = solve_daily_opf(data, p, q, pv, base_fixed)
        optimized = solve_daily_opf(data, p, q, pv)
        improvement = np.maximum(np.asarray(baseline["objective_t"]) - np.asarray(optimized["objective_t"]), 0.0)
        weights = improvement / improvement.sum() if improvement.sum() > 1e-12 else np.ones(data.n_time) / data.n_time
        action = np.abs(np.asarray(optimized["z"]) - BASE_STATUS[:, None])
        # Gao Eq. (19)--(23) analogue: performance share times action indicator,
        # accumulated over time and averaged over M operating samples.
        scores += action @ weights
        trace.append({"sample": sample_id, "baseline_objective": baseline["objective"],
                      "optimized_objective": optimized["objective"], "hourly_weight": weights.tolist()})
    return scores / sample_count, trace


def fcm_period_division(case_dir: Path, cluster_count: int, time_feature_weight: float,
                        seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Implement Gao et al. (2022) FCM time clustering (Eqs. 47--51).

    FCM is applied to the 16-dimensional hourly load vector.  A time feature
    ``time_feature_weight * t`` is appended, matching the existing FCM2.py
    implementation in this project and making the resulting operating periods
    consecutive, as required by the cited paper.  This continuity augmentation
    is reported explicitly; it is an implementation choice, not an equation in
    Gao et al. (2022).
    """
    raw_load = pd.read_excel(case_dir / "load_P.xlsx").values
    features = np.column_stack((np.arange(raw_load.shape[1]) * time_feature_weight, raw_load.T))
    np.random.seed(seed)
    _, memberships, _, _, _, _, _ = fuzz.cluster.cmeans(
        features.T, cluster_count, 2.0, error=1e-6, maxiter=1000, init=None
    )
    return np.argmax(memberships, axis=0), memberships


def write_outputs(case_dir: Path, output_dir: Path, sample_count: int, key_count: int, seed: int,
                  period_count: int, time_feature_weight: float) -> None:
    data = read_case(case_dir)
    if data.n_branch != len(BASE_STATUS):
        raise ValueError("This reference implementation is defined for the 16-branch IEEE three-feeder case.")
    scores, trace = contribution_scores(data, sample_count, seed)
    order = np.argsort(-scores, kind="stable")
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "switch_contributions.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "branch_index_0_based", "from_bus", "to_bus", "base_status", "contribution_score"])
        for rank, branch in enumerate(order, start=1):
            writer.writerow([rank, branch, data.branch_in[branch], data.branch_out[branch], BASE_STATUS[branch], f"{scores[branch]:.10f}"])
    with (output_dir / "key_switches.txt").open("w", encoding="utf-8") as f:
        f.write("Static key switches (zero-based branch indices): " + ", ".join(map(str, order[:key_count])) + "\n")
        f.write("Static key switches (one-based branch indices): " + ", ".join(map(str, order[:key_count] + 1)) + "\n")
        f.write("Gao-style static key-switch set; FCM time division is reported separately.\n")
    with (output_dir / "run_metadata.txt").open("w", encoding="utf-8") as f:
        f.write(f"samples={sample_count}\nseed={seed}\nkey_count={key_count}\n")
        f.write(f"fcm_period_count={period_count}\nfcm_time_feature_weight={time_feature_weight}\n")
        for row in trace:
            f.write(f"sample={row['sample']}, baseline={row['baseline_objective']:.8f}, optimized={row['optimized_objective']:.8f}\n")
    labels, memberships = fcm_period_division(case_dir, period_count, time_feature_weight, seed)
    period_ids = np.zeros(data.n_time, dtype=int)
    period_ids[0] = 1
    for t in range(1, data.n_time):
        period_ids[t] = period_ids[t - 1] + int(labels[t] != labels[t - 1])
    with (output_dir / "fcm_time_periods.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["hour_0_based", "fcm_cluster", "contiguous_period", *[f"membership_cluster_{c}" for c in range(period_count)]])
        for t in range(data.n_time):
            writer.writerow([t, labels[t], period_ids[t], *[f"{memberships[c, t]:.10f}" for c in range(period_count)]])
    periods: list[tuple[int, int]] = []
    start = 0
    for t in range(1, data.n_time + 1):
        if t == data.n_time or labels[t] != labels[t - 1]:
            periods.append((start, t - 1))
            start = t
    with (output_dir / "fcm_time_periods.txt").open("w", encoding="utf-8") as f:
        f.write("FCM contiguous operating periods (zero-based hours): " + "; ".join(f"{a}-{b}" for a, b in periods) + "\n")
        f.write("FCM contiguous operating periods (clock time): " + "; ".join(f"{a:02d}:00-{b + 1:02d}:00" for a, b in periods) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Gao-style static switch-contribution baseline")
    parser.add_argument("--case-dir", type=Path, required=True, help="Directory containing branch.xlsx, load_P.xlsx, load_Q.xlsx and PV.xlsx")
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--samples", type=int, default=20, help="Number of reproducible operating samples, including the nominal trace")
    parser.add_argument("--key-count", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20260821)
    parser.add_argument("--period-count", type=int, default=5, help="FCM cluster count (five corresponds to at most four daily topology changes)")
    parser.add_argument("--fcm-time-weight", type=float, default=200.0, help="Weight of the appended hour coordinate used to produce contiguous operating periods")
    args = parser.parse_args()
    if args.samples < 1 or not 1 <= args.key_count <= 16:
        raise ValueError("samples must be >= 1 and key-count must be in 1..16")
    write_outputs(args.case_dir, args.output_dir, args.samples, args.key_count, args.seed,
                  args.period_count, args.fcm_time_weight)


if __name__ == "__main__":
    main()

# Gao-style static switch-contribution baseline

This folder is a reproducible, model-based comparison baseline inspired by the
switch-contribution screening in:

H. Gao, S. Jiang, Z. Li, *et al.*, “A Two-Stage Multi-Agent Deep
Reinforcement Learning Method for Urban Distribution Network Reconfiguration
Considering Switch Contribution,” *IEEE Transactions on Power Systems*, vol.
39, no. 6, pp. 7064–7076, 2024.

## What is reproduced and what is adapted

The reference first collects optimized switch actions under multiple operating
samples, then accumulates a performance weight whenever a switch changes state
(its Eqs. (17)–(23)), and finally retains the high-ranked switches as a static
set.  The original article uses PV-curtailment and load-shedding capability as
the performance weight and uses WQMIX to learn the subsequent reconfiguration
policy.

This implementation retains the static contribution-screening mechanism but
uses the objectives of the present study: network loss and branch-current load
imbalance.  For operating sample m and hour t, the performance weight is

`q(m,t) = max(J_base(m,t) - J_full(m,t), 0) / sum_t max(J_base(m,t) - J_full(m,t), 0)`.

`J_base` is the objective with the original 13 sectionalizing switches closed
and 3 tie switches open; `J_full` is the full-switch OPF objective.  The static
contribution of branch b is

`E_b = (1/M) sum_m sum_t q(m,t) |z_b(m,t) - z_b_base|`.

This is the direct analogue of Gao et al.’s action-weighted contribution, with
an objective-aligned weight.  It is not a claim to reproduce their DRL policy,
nor is it a period-division method.  Consequently, the result is **one static
key-switch set and no time partition**.

## Reproduction

Run from this folder with an environment containing `gurobipy`, `numpy`, and
`pandas`:

```powershell
python gao_static_contribution.py `
  --case-dir "..\5.IEEE典型三馈线网络 - 得到分段关键开关后的代码\trans_case\ieee16_1" `
  --output-dir results --samples 20 --key-count 6 --seed 20260821
```

The nominal daily trace is always sample 0.  The remaining 19 samples apply
small, seeded daily load/PV variations, so results are deterministic.  The
program writes `switch_contributions.csv`, `key_switches.txt`,
`fcm_time_periods.csv`, `fcm_time_periods.txt`, and `run_metadata.txt` to the
selected output directory.

## FCM period division

The FCM portion follows Gao et al., “Multi-objective Dynamic Reconfiguration
for Urban Distribution Network Considering Multi-level Switching Modes,”
*Journal of Modern Power Systems and Clean Energy*, 10(5), 1241–1255, 2022.
It uses the FCM objective and membership updates in its Eqs. (47)–(50), with
fuzziness exponent 2 and five clusters. Five clusters are consistent with the
present study’s limit of at most four topology changes per day.

The article requires every cluster to correspond to a continuous operating
period. To obtain that property on this 24-hour data set, the code explicitly
appends the hour coordinate `200*t` to the raw spatial load vector; this is the
same continuity device already used in this project’s `FCM2.py`. It is an
implementation choice disclosed in `run_metadata.txt`, not an unstated claim
about Gao et al.’s equations. The resulting time periods and the key switches
are deliberately decoupled and static, which is the intended baseline.

## Fair comparison statement

Use the same 16-node data, the same 24-hour horizon, the same radiality and
power-flow constraints, and the same loss:imbalance objective ratio (1:7) as
the proposed method.  The comparison isolates the distinction of interest:
Gao-style **static spatial selection + FCM period division** versus the
proposed **period-adaptive spatiotemporal selection**. The cited 2024 switch
contribution paper itself does not divide time periods; the FCM component is
instead attributed to Gao et al. (2022), as stated above.

## CP-SAT for fixed-$K$ incumbent generation

CP was considered a potentially faster alternative to MILP for the fixed-$K$ feasibility check. We therefore implemented CP-SAT backends that construct exactly $K$ routes and convert every accepted solution into a MIP start for the existing multigraph MILP. Three variants were compared.

- **OG-S**: an anchored-circuit **satisfaction model** in which pickup, treatment, and delivery are represented as atomic actions. It has no objective and directly enforces $C_k\le T$ while searching for the first feasible solution.
- **OG-M($\rho$)**: the same atomic-action model with the relaxed horizon $U=\lceil\rho T\rceil$ and **objective $\min\max_k C_k$**. An incumbent satisfying $\max_k C_k\le T$ is a feasible warm start for the original problem.
- **MG-M($\rho$)**: a **multigraph-based CP** model that uses the existing state-expanded multigraph edges as circuit arcs and applies the **same makespan objective.** Treatment operations and skip-state transitions are embedded in the edges, but the edge-specific Boolean literals and channeling constraints require substantial memory.
- **+R**: the **corresponding redundant-strengthening bundle**, including terminal balance, full-skip tracking, container workload, and aggregate-duration constraints.

Each exact-$K$ attempt has a 1200-second time limit. All entries are wall-clock seconds. `U` denotes `UNKNOWN` without an incumbent, `I` denotes proven infeasibility, and `KILL` denotes an out-of-memory termination by operating-system signal 9.

| ID | $K$ | MILP duration | OG-S | OG-S+R | OG-M(1.5) | OG-M(1.5)+R | MG-M(1.5)+R |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| A1 | 1 | 0.02 | 0.05 | 0.06 | 0.05 | 0.05 | 0.14 |
| A2 | 1 | 0.02 | 0.09 | 0.13 | 0.11 | 0.13 | 0.34 |
| A3 | 1 | 0.05 | 0.11 | 0.12 | 0.09 | 0.13 | 0.73 |
| A4 | 1 | 0.07 | 136.82 | 188.93 | 2.26 | 2.81 | 2.38 |
| A5 | 1 | 0.11 | U (1039.19) | U (221.43) | 3.69 | 4.03 | 3.20 |
| A6 | 2 | 0.04 | 10.57 | 4.23 | 4.54 | 2.86 | 8.79 |
| A7 | 2 | 0.13 | 7.83 | 5.71 | 2.45 | 2.65 | 7.28 |
| **A1--A7 total** |  | **0.44** | **1194.66** | **420.61** | **13.19** | **12.66** | **22.86** |
| D5 | 7 | I (126.71) | -- | -- | U (1220.39) | -- | KILL |

## Acceleration of the duration IP and incumbent search

Two acceleration mechanisms were introduced for the duration-based lower-bound problem. **First, among parallel multigraph edges having the same tail node, head node, and endpoint states, only a minimum-duration representative is retained**. Second, **internal delivery-to-pickup connectors whose intermediate state is empty are removed:** because the duration IP does not fix the number of routes, such a connector can be split into two open routes without increasing total duration. The second elimination is not valid for the exact-$K$ incumbent problem and is therefore used only for the duration IP.

**The duration IP is also terminated as soon as its incumbent and global bound imply the same rounded fleet bound,**

$$
\left\lceil\frac{H^{\mathrm{LB}}}{T}\right\rceil
=
\left\lceil\frac{H^{\mathrm{UB}}}{T}\right\rceil.
$$

This certifies the exact value of $\lceil H^*/T\rceil$ without proving the optimal duration $H^*$ itself. The fixed-$K$ incumbent model uses only minimum-duration parallel-edge elimination and still terminates at the first validated feasible solution.

### Duration-IP results

Every accelerated run terminated through the rounded-bound callback and reproduced the previously obtained duration-based fleet bound $\widehat k_D$. `RB` denotes rounded-bound certification, and `TL` denotes the original 600-second time limit.

| Instance | $\widehat k_D$ | Original time/status | Accelerated time/status | Recorded speedup |
| :-- | --: | --: | --: | --: |
| A1--A7 | 1--2 | 1.06 s total | 0.43 s total / RB | 2.5$\times$ |
| C16 | 8 | 600.11 s / TL | 12.12 s / RB | $\ge49.5\times$ |
| C17 | 9 | 43.91 s | 14.75 s / RB | 3.0$\times$ |
| C18 | 9 | 600.28 s / TL | 15.90 s / RB | $\ge37.7\times$ |
| C19 | 12 | 600.01 s / TL | 18.04 s / RB | $\ge33.3\times$ |
| C20 | 12 | 600.04 s / TL | 25.26 s / RB | $\ge23.8\times$ |
| D3 | 5 | 30.32 s | 3.89 s / RB | 7.8$\times$ |
| D4 | 7 | 15.19 s | 5.98 s / RB | 2.5$\times$ |
| D5 | 7 | 600.02 s / TL | 10.43 s / RB | $\ge57.5\times$ |
| D6 | 7 | 21.09 s | 9.38 s / RB | 2.2$\times$ |
| D7 | 9 | 51.83 s | 12.30 s / RB | 4.2$\times$ |

Across these 17 instances, recorded duration-IP time decreased from 3163.86 to 128.49 seconds, a factor of 24.6. The accelerated duration graphs contained 2.91--4.30% fewer edges than the corresponding main graphs. Since arc elimination and rounded-bound termination were enabled together, this table measures their combined effect.

### Fixed-$K$ parallel-edge elimination

The following controlled comparison uses the revised treatment-state transition and a 1200-second limit in both runs. T**he incumbent graph with arc elimination removes only dominated parallel edges**; empty-state connector elimination remains disabled. **Across all 67 instances, 0.81% of incumbent edges were removed in aggregate, with an instance-level range of 0--2.22%.**

The next table aggregates groups for which both configurations produced the same feasibility outcome. `F` denotes a validated incumbent and `U` denotes time-limit termination without an incumbent.

| Instances | Common outcome | Without elimination | With elimination |
| :-- | :--: | --: | --: |
| A1--A20 | 20 F | 56.93 s total | 73.22 s total |
| B1--B20 | 20 F | 115.36 s total | 140.49 s total |
| C1--C11 | 11 F | 509.56 s total | 321.04 s total |
| C12, C16, C18, C20 | 4 U | 4800.37 s total | 4800.36 s total |
| C13--C15 | 3 F | 341.09 s total | 262.70 s total |
| D1--D5 | 5 F | 322.47 s total | 294.98 s total |
| D7 | 1 F | 233.77 s | 140.71 s |

The feasibility outcome differed on only three instances. The reported UB is the duration of the first incumbent, while LB is the final global duration bound.

| Instance | $K$ |     Edges: without $\rightarrow$ with | Without elimination            | With elimination                |
| :------- | --: | ------------------------------------: | :----------------------------- | :------------------------------ |
| C17      |   9 | $231071\rightarrow229580$ ($-0.65\%$) | U / 1200.10 s, LB 4185         | F / 1098.08 s, UB 4202, LB 4183 |
| **C19**  |  12 | $363699\rightarrow361312$ ($-0.66\%$) | F / 476.13 s, UB 5559, LB 5475 | U / 1200.27 s, LB 5477          |
| D6       |   7 | $230630\rightarrow228033$ ($-1.13\%$) | U / 1200.04 s, LB 3222         | F / 1022.83 s, UB 3308, LB 3219 |

Overall, the outcome count changed from 61 F and 6 U without elimination to 62 F and 5 U with elimination. Among the 60 instances solved feasibly by both configurations, total runtime decreased from 1579.18 to 1233.14 seconds, but the effect was not uniform: elimination was faster on 26 instances, slower on 24, and tied on 10. In particular, it enabled incumbents for C17 and D6 but lost the incumbent previously found for C19. The experiment therefore supports parallel-edge elimination as a safe model reduction with occasional large search benefits, but not as a monotone runtime improvement.

## Revision of the treatment-location state transition

### Modeling issue and correction

**The original transition emptied every onboard full skip assigned to a visited treatment location.** When two eligible full skips were onboard, it therefore excluded the valid decision to empty only one skip and leave the other full until a later visit. Treating the all-at-once decision as dominant implicitly requires a metric travel-time matrix. **The current instances do not satisfy the triangle inequality**, **so removing the later treatment visit can increase route duration and can remove feasible or lower-cost routes.**


For example, suppose two skips use the same treatment location $h$ and a feasible action sequence contains

$$
P_1\rightarrow P_2\rightarrow H_1(h)\rightarrow D_1
\rightarrow H_2(h)\rightarrow D_2.
$$

If both skips are forced to be emptied at the first visit to $h$, the second visit is removed and the corresponding segment becomes $D_1\rightarrow D_2$. With

$$
t_{D_1,h}=1,\qquad t_{h,D_2}=1,\qquad t_{D_1,D_2}=10,
$$

this replacement increases travel time from $2$ to $10$. Under the triangle inequality this cannot occur because $t_{D_1,D_2}\le t_{D_1,h}+t_{h,D_2}$. The revised partial-emptying transition is therefore necessary for the non-metric instances used in this project.

The revised transition retains the existing state representation but **generates three outcomes whenever two full skips are eligible at the visited treatment: empty the first skip, empty the second skip, or empty both.** If only one skip is eligible, it is emptied as before. The correction changes only the edge set by adding parallel edges for the partial-emptying outcomes.
### Invalidity of the original COR workload bound

The same non-metric travel times invalidate the original COR workload expression

$$
k_{\min}^{\mathrm{COR}}
=
\left\lceil
\frac{
\left(
\displaystyle\sum_{r\in R} t_{\alpha(r),\beta(r)}
+
\displaystyle\sum_{r\in R}
\min_{\substack{r'\in R\\ \phi(r')=\phi(r)}}
t_{\beta(r),\omega(r')}
\right)/2
+
\displaystyle\sum_{i\in N}s_i
}{T}
\right\rceil.
$$

Its travel component assumes that the direct entry $t_{uv}$ is no greater than the travel time of any path from $u$ to $v$. Without the triangle inequality, the direct pickup-to-treatment and treatment-to-delivery entries need not be lower bounds on the corresponding travel embedded in an interleaved route.

Consider two same-type requests and the route order

$$
\alpha_1\rightarrow\alpha_2\rightarrow\beta_1
\rightarrow\beta_2\rightarrow\omega_1\rightarrow\omega_2,
$$

where each consecutive arc has travel time $1$. The route travel time is therefore $5$. However, let

$$
t_{\alpha_1,\beta_1}
=t_{\alpha_2,\beta_2}
=t_{\beta_1,\omega_1}
=t_{\beta_2,\omega_2}
=100.
$$

The COR travel term then becomes

$$
\frac{100+100+100+100}{2}=200,
$$

despite the existence of a route with travel time $5$. Hence the expression can overestimate the required workload and is **not a valid vehicle lower bound for the current non-metric instances**. The previously computed COR values must not be treated as certified lower bounds unless the formula is rederived using valid metric-closure distances or replaced by another proven bound.

### Graph-size effect

The state count is unchanged in all 67 instances. The final edge count increases by **2.70% on average**, with a range of **0--7.64%** across instances.

### Duration-based fleet bound

The selected duration-based fleet bound $\widehat k_D$ is unchanged on all 67 instances. 

### Fixed-$K$ incumbent search

The existing runs use a 600-second limit, whereas the revised runs use 1200 seconds and the specialized incumbent graph. Consequently, this is a descriptive before-and-after comparison rather than an isolated state-transition ablation. `F`, `I`, and `U` denote a validated feasible incumbent, proven infeasibility, and unresolved time-limit termination, respectively.

The following groups have the same outcome under both formulations.

| Instances | Common outcome | Existing time | Revised time |
| :-- | :--: | --: | --: |
| A1--A20 | 20 F | 60.69 s total | 73.22 s total |
| B1--B20 | 20 F | 120.04 s total | 140.49 s total |
| C1--C11 | 11 F | 360.01 s total | 321.04 s total |
| C13--C15 | 3 F | 241.13 s total | 262.70 s total |
| C16, C18, C20 | 3 U | 1800.44 s total | 3600.34 s total |
| D1--D4, D7 | 5 F | 167.52 s total | 183.71 s total |

The outcome changes are as follows. The reported cost is the original-cost value of the accepted fixed-$K$ incumbent.

| Instance | $K$ | Existing formulation | Revised formulation |
| :-- | --: | :-- | :-- |
| C12 | 5 | F / 288.24 s, cost 3849 | U / 1200.02 s |
| C17 | 9 | U / 600.07 s | F / 1098.08 s, cost 7066 |
| C19 | 12 | I / 571.36 s | U / 1200.27 s |
| D5 | 7 | I / 126.71 s | F / 251.98 s, cost 5047 |
| D6 | 7 | U / 600.17 s | F / 1022.83 s, cost 5018 |

### Main 2-index IP

The main-IP comparison uses a common 600-second cutoff for the smaller set and a common 300-second cutoff for the larger set. At these cutoffs, 46 of 67 instances have exactly the same UB, LB, and gap. These are

| Identical results | Number of instances |
| :-- | --: |
| A1--A12, A14--A16, A18--A20; B2--B5, B7--B20; C1--C5, C7, C10; D1, D2, D4 | 46 |

The remaining instances are reported individually below.

| Instance | Existing UB / LB / gap (%) | Revised UB / LB / gap (%) |
| :------- | :------------------------- | :------------------------ |
| A13      | 1379 / 1375 / 0.29         | 1378 / 1375 / 0.22        |
| **A17**  | **2101 / 2101 / 0.00**     | **2112 / 2099 / 0.62**    |
| B1       | 1291 / 1284.48 / 0.50      | 1291 / 1285 / 0.46        |
| B6       | 1171 / 1166 / 0.43         | 1172 / 1166 / 0.51        |
| C6       | 3009 / 2987 / 0.73         | 3008 / 2986 / 0.73        |
| C8       | 3044 / 3006 / 1.25         | 3032 / 3010 / 0.73        |
| C9       | 3960 / 3916.44 / 1.10      | 3931 / 3917 / 0.36        |
| C11      | 4650 / 4650 / 0.00         | 4632 / 4632 / 0.00        |
| C12      | 3844 / 3792.50 / 1.34      | 3800 / 3799 / 0.03        |
| C13      | 5313 / 5257 / 1.05         | 5386 / 5256.38 / 2.41     |
| **C14**  | **5462 / 5462 / 0.00**     | **5486 / 5453.67 / 0.59** |
| C15      | 6277 / 6166 / 1.77         | 6242 / 6162 / 1.28        |
| C16      | 6865 / 6235 / 9.18         | 6817 / 6242 / 8.43        |
| C17      | 7571 / 7026 / 7.20         | 7059 / 7042.50 / 0.23     |
| C18      | 7694 / 7045.01 / 8.44      | 7621 / 7042.56 / 7.59     |
| C19      | 10823 / 9367 / 13.45       | 10908 / 9366.82 / 14.13   |
| C20      | 9987 / 9366 / 6.22         | 10122 / 9367.81 / 7.45    |
| **D3**       | **3522 / 3514 / 0.23**         | **3511 / 3511 / 0.00**        |
| D5       | 5534 / 4910 / 11.28        | 4998 / 4908.78 / 1.79     |
| D6       | 5391 / 4885 / 9.39         | 4995 / 4886 / 2.18        |
| D7       | 6425 / 6292 / 2.07         | 6350 / 6292.82 / 0.90     |

The revised model loses an optimality proof within the cutoff for A17 and C14, but proves D3 optimal. More importantly, C11 obtains a strictly better proven optimum, decreasing from 4650 to 4632. This confirms that the original transition could exclude valid lower-cost solutions, whereas the revised transition represents the intended SPDP treatment decisions.

## DFF-based route-bound enhancement

### Two uses of dual-feasible functions

A dual-feasible function (DFF) $f:[0,1]\rightarrow[0,1]$ satisfies

$$
\sum_i x_i\le1
\quad\Longrightarrow\quad
\sum_i f(x_i)\le1.
$$

Since the normalized edge durations on every feasible route satisfy $\sum_{e\in r}t_e/T\le1$, a DFF yields the valid route inequality

$$
\boxed{
K(\theta)
\ge
\sum_{e\in E}
f\!\left(\frac{t_e}{T}\right)\theta_e
}.
$$

This validity follows directly from route-duration feasibility. DFFs were tested in two different roles.

1. **DFF-weighted Sub-IPs:** replace the duration objective by
   $$
   \min\sum_{e\in E}f(t_e/T)y_e
   $$
   and use the ceiling of its certified bound as a candidate fixed route lower bound.
2. **Direct incumbent cuts:** add the boxed inequality to the fixed-$K$ duration MILP used to generate the warm start.

The tests used the identity function and the Fekete--Schepers family with $\lambda\in\{0.1,0.2,0.3,0.4\}$. The identity case reproduces the ordinary normalized duration bound. The strongest Fekete--Schepers result was always obtained with $\lambda=0.1$.

### DFF-weighted duration Sub-IPs

The DFF-weighted Sub-IPs did not improve the route lower bound on any instance. Among the 66 completed runs, the strongest Fekete--Schepers bound was equal to the ordinary duration bound on only three small instances and was strictly weaker on the other 64.

| Comparison with the ordinary duration bound | Instances                                 | Count |
| :------------------------------------------ | :---------------------------------------- | ----: |
| Equal                                       | A1, A3, A5                                |     3 |
| Strictly weaker                             | A2, A4, A6--A20; B1--B20; C1--C20; D1--D7 |    64 |
| Stronger                                    | --                                        |     0 |


Consequently, the final $k_{\min}$ remained unchanged on all 67 instances because the ordinary duration IP dominated every completed nonlinear DFF result. The additional DFF Sub-IPs required approximately 1744.74 seconds, including 600-second time-limit runs for $\lambda=0.1$ on C9 and C19. In this edge-duration application, the tested DFFs therefore produced bounds that were only equal or worse while adding substantial computational cost.

### Direct DFF cuts in the fixed-$K$ incumbent MILP

The incumbent experiment added five inequalities: one identity cut and four Fekete--Schepers cuts. Each entry below is based on SolutionLimit = 1; the UB is the duration of the first validated incumbent and the LB is the final global bound. F, TL, and OOM denote feasible, time-limit without an incumbent, and out-of-memory termination, respectively.

The following groups had the same outcome with and without the direct cuts.

| Instances | Common outcome | Without cuts | With cuts |
| :-- | :--: | --: | --: |
| A1--A20 | 20 F | 73.22 s total | 79.91 s total |
| B1--B20 | 20 F | 140.49 s total | 112.72 s total |
| C1--C11 | 11 F | 321.04 s total | 178.52 s total |
| C13--C15 | 3 F | 262.70 s total | 301.32 s total |
| C18, C20 | 2 TL | 2400.28 s total | 2400.28 s total |
| D1--D5, D7 | 6 F | 435.69 s total | 97.91 s total |

**The outcome changed on five instances.**

| Instance | $K$ | Without direct DFF cuts | With direct DFF cuts |
| :-- | --: | :-- | :-- |
| C12 | 5 | TL / 1200.02 s, LB 2332 | F / 637.87 s, UB 2380, LB 2328 |
| C16 | 8 | TL / 1200.06 s, LB 3705 | F / 422.33 s, UB 3727, LB 3707 |
| C17 | 9 | F / 1098.08 s, UB 4202, LB 4183 | TL / 1200.12 s, LB 4174 |
| C19 | 12 | TL / 1200.27 s, LB 5477 | F / 69.45 s, UB 5523, LB 5450 |
| D6 | 7 | F / 1022.83 s, UB 3308, LB 3219 | OOM |

Without the cuts, the complete 67-instance experiment produced 62 feasible incumbents and five time-limit outcomes. With the cuts, it produced 63 feasible incumbents, three time-limit outcomes, and one OOM outcome. On the 66-instance paired subset excluding D6, total incumbent time decreased from 8331.85 to 5500.43 seconds, a reduction of approximately 34.0%.

The direct inequalities were therefore useful as search-strengthening constraints, particularly for C12, C16, and C19, but their effect was not uniform: C17 changed from feasible to time-limit and D6 exhausted memory. The two DFF applications should consequently be judged separately. The DFF-weighted Sub-IPs are not supported by the present evidence, whereas the direct incumbent cuts remain potentially useful as an optional strengthening.

## Paper: Porumbel, D., & Goncalves, G. (DAM, 2015). Using dual feasible functions to construct fast lower bounds for routing and location problems (DFF-based dual bounds for arc-routing variants)

### Main idea

Porumbel and Goncalves (2015) construct fast lower bounds for route-column formulations by generating a **feasible solution of the master dual directly from dual-feasible functions**, rather than solving the complete root-node LP by column generation. For a route $r$, a representative dual constraint has the form

$$
\sum_{j\in E} a_{jr}y_j-\mu\le c_r,
\qquad r\in\mathcal R,
$$

where $a_{jr}$ is the number of services of edge $j$ in route $r$. **A DFF maps normalized capacity or distance consumption to dual weights while guaranteeing that their sum is at most one on every feasible route.** The resulting dual variables are restricted to a low-dimensional parameterized family, so only a small auxiliary optimization problem must be solved and all exponentially many route constraints are satisfied analytically.

This gives the following hierarchy.

| Method | Treatment of route dual constraints | Bound quality | Computational effort |
| :-- | :-- | :-- | :-- |
| Pure DFF | Guaranteed by the DFF construction | Weakest | Very small |
| Mixed CG--DFF | Selected edge duals remain unrestricted | Intermediate | Intermediate |
| Full column generation | Pricing enforces the complete route dual | Exact root-node LP bound | Largest |

Thus, the DFF bound cannot dominate the exact column-generation root bound, but it can provide a reasonably strong lower bound almost immediately. It can also be used as a warm start for column generation.

### Deadheading and bound loss

Let $a_{jr}^{\mathrm{tr}}$ denote the total number of traversals of edge $j$, including deadheading. Since route cost contains every traversal whereas the covering dual rewards serviced edges, substitution of the DFF dual values leaves the nonnegative slack

$$
\sum_{j\in E}
\left(a_{jr}^{\mathrm{tr}}-a_{jr}\right)d_j.
$$

**The pure DFF construction discards this route-dependent slack to obtain a simple sufficient condition that is valid for every route.** Consequently, it cannot fully price the cost of reaching required edges or traveling between services. Its quality deteriorates when deadheading forms a large fraction of route cost, as observed particularly on large CARP instances containing many non-required edges.

### Why distance-constrained ARP is more favorable

In the distance-constrained ARP, route feasibility is defined directly by total traversal distance,

$$
\sum_{j\in E}a_{jr}^{\mathrm{tr}}d_j\le D.
$$

Moreover, the route-column master uses covering constraints and does not impose an upper bound on the number of times an edge may be serviced. Therefore, any traversal currently labeled as deadheading can be relabeled as service without changing the route, its cost, or its distance feasibility. For the dual constraints that can become tight, one may consequently assume

$$
a_{jr}=a_{jr}^{\mathrm{tr}}
\qquad \forall j\in E.
$$

**The deadheading slack then vanishes, and the distance DFF is applied to the actual traversal workload, rather than only to the serviced portion of a route.** This alignment between route cost, resource consumption, and dual incidence explains why DFF bounds are substantially stronger for the distance-constrained ARP than for general capacitated arc routing.

On the standard distance-constrained instances with fewer than 50 edges, the pure DFF bound reached approximately **90--99% of the column-generation bound while requiring only about 0.1--0.5% of its computation time**. On the larger `egl` instances, the bound remained very fast but decreased to roughly 75% of the CG bound. The principal contribution is therefore not an exact replacement for column generation, but a fast, analytically dual-feasible bound whose quality is particularly good when deadheading can be absorbed into the resource representation.

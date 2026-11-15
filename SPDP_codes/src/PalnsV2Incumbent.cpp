#include "PalnsV2Incumbent.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <thread>
#include <tuple>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "PstepFormulation.h"
#include "RouteState.h"

namespace spdp {

namespace {

using Clock = std::chrono::steady_clock;
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kEps = 1e-7;
constexpr std::int8_t kPick = 0;
constexpr std::int8_t kEmpty = 1;
constexpr std::int8_t kDeliv = 2;

struct Op {
    std::int8_t kind = kPick;
    int cls = -1;   // pickup class for P/E, delivery class for D
    int traj = -1;  // skip trajectory id
};

struct Instance {
    int n = 0;
    double sp = 0.0, se = 0.0, sd = 0.0, T = 0.0, F = 0.0;
    const std::vector<std::vector<double>>* t = nullptr;
    const std::vector<std::vector<double>>* c = nullptr;
    std::vector<int> pa, pb, pf, pm;
    std::vector<std::vector<int>> pmem;
    std::vector<int> dw, df, dm;
    std::vector<std::vector<int>> dmem;
    std::vector<int> req_pc, req_dc;
    int ntype = 0;

    double tt(int a, int b) const { return (*t)[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]; }
    double cc(int a, int b) const { return (*c)[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]; }
    double lt(int a, int b) const { return (a < 0 || b < 0) ? 0.0 : tt(a, b); }
    double lc(int a, int b) const { return (a < 0 || b < 0) ? 0.0 : cc(a, b); }
    int loc(const Op& o) const {
        return o.kind == kPick ? pa[static_cast<std::size_t>(o.cls)]
             : (o.kind == kEmpty ? pb[static_cast<std::size_t>(o.cls)]
                                 : dw[static_cast<std::size_t>(o.cls)]);
    }
    double srv(const Op& o) const {
        return o.kind == kPick ? sp : (o.kind == kEmpty ? se : sd);
    }
    double single_d(int p, int q) const {
        return sp + se + sd + tt(pa[static_cast<std::size_t>(p)], pb[static_cast<std::size_t>(p)]) +
            tt(pb[static_cast<std::size_t>(p)], dw[static_cast<std::size_t>(q)]);
    }
    double single_c(int p, int q) const {
        return cc(pa[static_cast<std::size_t>(p)], pb[static_cast<std::size_t>(p)]) +
            cc(pb[static_cast<std::size_t>(p)], dw[static_cast<std::size_t>(q)]);
    }
};

inline int load_delta(const Op& o) {
    return o.kind == kPick ? 1 : (o.kind == kDeliv ? -1 : 0);
}

inline double excess(double d, double T) {
    return d > T ? d - T : 0.0;
}

struct Route {
    std::vector<Op> ops;
    double d = 0.0;
    double c = 0.0;
};

struct Solution {
    std::vector<Route> routes;
    double V = 0.0;
    double D = 0.0;
    double C = 0.0;
};

void eval_ops(const Instance& I, const std::vector<Op>& ops, double& d, double& c) {
    d = 0.0;
    c = 0.0;
    int prev = -1;
    for (const Op& o : ops) {
        const int l = I.loc(o);
        d += I.srv(o);
        if (prev >= 0) {
            d += I.tt(prev, l);
            c += I.cc(prev, l);
        }
        prev = l;
    }
}

void refresh_route(const Instance& I, Route& r) {
    eval_ops(I, r.ops, r.d, r.c);
}

void refresh_totals(const Instance& I, Solution& s) {
    s.V = 0.0;
    s.D = 0.0;
    s.C = 0.0;
    for (const Route& r : s.routes) {
        s.V += excess(r.d, I.T);
        s.D += r.d;
        s.C += r.c;
    }
}

int nonempty_route_count(const Solution& s) {
    int count = 0;
    for (const Route& r : s.routes) {
        if (!r.ops.empty()) {
            ++count;
        }
    }
    return count;
}

void route_loads(const std::vector<Op>& ops, std::vector<int>& L) {
    L.assign(ops.size() + 1U, 0);
    for (std::size_t i = 0; i < ops.size(); ++i) {
        L[i + 1U] = L[i] + load_delta(ops[i]);
    }
}

// ---------------------------------------------------------------- patterns

struct SingleOpt {
    int q = -1;
    double d = 0.0;
    double c = 0.0;
};

struct PatOp {
    std::int8_t kind = kPick;
    std::int8_t slot = 0;  // P/E: pickup slot, D: delivery slot in order
};

struct Composite {
    int p0 = -1, p1 = -1, q0 = -1, q1 = -1;  // q0: first D, q1: second D
    std::int8_t skel = 0;
    std::array<PatOp, 6> ops{};
    std::array<std::int8_t, 2> downer{{0, 1}};  // pickup slot owning each D
    double d = 0.0, c = 0.0, delta = 0.0;
};

struct PatternPool {
    std::vector<std::vector<SingleOpt>> singles;
    std::vector<Composite> comps;
    std::unordered_map<long long, std::vector<int>> comp_by_pair;
    std::vector<std::vector<double>> sim;
    std::array<double, 4> sim_raw_range{{0.0, 0.0, 0.0, 0.0}};  // lo/hi of the two raw travel terms
    std::size_t raw = 0, duplicates = 0, dominated = 0;
};

long long pair_key(int a, int b, int P) {
    if (a > b) {
        std::swap(a, b);
    }
    return static_cast<long long>(a) * static_cast<long long>(P) + b;
}

void build_pattern_pool(const Instance& I, const PalnsOptions& opt, PatternPool& pool) {
    const int P = static_cast<int>(I.pa.size());
    const int Q = static_cast<int>(I.dw.size());
    pool.singles.assign(static_cast<std::size_t>(P), {});
    for (int p = 0; p < P; ++p) {
        for (int q = 0; q < Q; ++q) {
            if (I.df[static_cast<std::size_t>(q)] != I.pf[static_cast<std::size_t>(p)]) {
                continue;
            }
            pool.singles[static_cast<std::size_t>(p)].push_back(
                SingleOpt{q, I.single_d(p, q), I.single_c(p, q)});
        }
        auto& list = pool.singles[static_cast<std::size_t>(p)];
        std::sort(list.begin(), list.end(), [](const SingleOpt& a, const SingleOpt& b) {
            if (a.d != b.d) return a.d < b.d;
            if (a.c != b.c) return a.c < b.c;
            return a.q < b.q;
        });
    }

    // SIM: the two travel-time terms are min-max normalized over all pairs
    // p != p' so that each term ranges over [0, 1] like the type indicator.
    const std::size_t PS = static_cast<std::size_t>(P);
    std::vector<std::vector<double>> raw1(PS, std::vector<double>(PS, 0.0)), raw2 = raw1;
    double lo1 = kInf, hi1 = -kInf, lo2 = kInf, hi2 = -kInf;
    for (std::size_t p = 0; p < PS; ++p) {
        for (std::size_t r = 0; r < PS; ++r) {
            raw1[p][r] = 1.0 / (1.0 + I.tt(I.pa[p], I.pa[r]));
            raw2[p][r] = 1.0 / (1.0 + I.tt(I.pb[p], I.pb[r]));
            if (p == r) continue;
            lo1 = std::min(lo1, raw1[p][r]); hi1 = std::max(hi1, raw1[p][r]);
            lo2 = std::min(lo2, raw2[p][r]); hi2 = std::max(hi2, raw2[p][r]);
        }
    }
    auto norm = [](double x, double lo, double hi) {
        if (!(hi > lo + 1e-12)) return x > hi + 1e-12 ? 1.0 : 0.0;
        return std::min(1.0, std::max(0.0, (x - lo) / (hi - lo)));
    };
    pool.sim.assign(PS, std::vector<double>(PS, 0.0));
    for (std::size_t p = 0; p < PS; ++p) {
        for (std::size_t r = 0; r < PS; ++r) {
            const double s1 = norm(raw1[p][r], lo1, hi1);
            const double s2 = norm(raw2[p][r], lo2, hi2);
            const double s3 = I.pf[p] == I.pf[r] ? 1.0 : 0.0;
            pool.sim[p][r] = opt.sim_s1 * s1 + opt.sim_s2 * s2 + opt.sim_s3 * s3;
        }
    }
    if (P > 1) {
        pool.sim_raw_range = {{lo1, hi1, lo2, hi2}};
    }

    const int kappa = std::max(1, opt.composite_kappa);
    auto lambda = [&](int p) {
        std::vector<int> out;
        for (const SingleOpt& s : pool.singles[static_cast<std::size_t>(p)]) {
            if (s.d <= I.T + kEps) {
                out.push_back(s.q);
                if (static_cast<int>(out.size()) >= kappa) {
                    break;
                }
            }
        }
        return out;
    };
    std::vector<std::vector<int>> lambdas(static_cast<std::size_t>(P));
    for (int p = 0; p < P; ++p) {
        lambdas[static_cast<std::size_t>(p)] = lambda(p);
    }
    const bool all_skeletons = opt.composite_skeletons == PalnsSkeletonSet::All;
    const double service2 = 2.0 * (I.sp + I.se + I.sd);

    for (int p0 = 0; p0 < P; ++p0) {
        for (int p1 = 0; p1 < P; ++p1) {
            if (p0 == p1 && I.pm[static_cast<std::size_t>(p0)] < 2) {
                continue;
            }
            const std::vector<int>& L0 = lambdas[static_cast<std::size_t>(p0)];
            const std::vector<int>& L1 = lambdas[static_cast<std::size_t>(p1)];
            std::set<std::pair<int, int>> dpairs;
            for (int qa : L0) {
                for (int qb : L1) {
                    if (qa == qb && I.dm[static_cast<std::size_t>(qa)] < 2) {
                        continue;
                    }
                    dpairs.insert(std::minmax(qa, qb));
                }
            }
            const int f0 = I.pf[static_cast<std::size_t>(p0)];
            const int f1 = I.pf[static_cast<std::size_t>(p1)];
            for (const auto& dp : dpairs) {
                const int qa = dp.first;
                const int qb = dp.second;
                std::vector<Composite> cands;
                std::set<std::vector<int>> seen;
                auto add = [&](std::int8_t skel, std::array<PatOp, 6> ops, int x, int y,
                               std::int8_t ox, std::int8_t oy) {
                    Composite cp;
                    cp.p0 = p0; cp.p1 = p1; cp.q0 = x; cp.q1 = y;
                    cp.skel = skel; cp.ops = ops; cp.downer = {{ox, oy}};
                    double travel = 0.0, cost = 0.0;
                    int prev = -1;
                    std::vector<int> signature;
                    for (const PatOp& po : ops) {
                        int l = -1;
                        int cls = -1;
                        if (po.kind == kPick) { cls = po.slot ? p1 : p0; l = I.pa[static_cast<std::size_t>(cls)]; }
                        else if (po.kind == kEmpty) { cls = po.slot ? p1 : p0; l = I.pb[static_cast<std::size_t>(cls)]; }
                        else { cls = po.slot ? y : x; l = I.dw[static_cast<std::size_t>(cls)]; }
                        if (prev >= 0) { travel += I.tt(prev, l); cost += I.cc(prev, l); }
                        prev = l;
                        signature.push_back(po.kind);
                        signature.push_back(l);
                        signature.push_back(po.kind == kDeliv ? cls : -1);
                    }
                    cp.d = service2 + travel;
                    cp.c = cost;
                    if (cp.d > I.T + kEps) {
                        return;
                    }
                    ++pool.raw;
                    if (!seen.insert(signature).second) {
                        ++pool.duplicates;
                        return;
                    }
                    cands.push_back(cp);
                };
                for (int ord = 0; ord < 2; ++ord) {
                    if (ord == 1 && qa == qb) {
                        break;
                    }
                    const int x = ord == 0 ? qa : qb;
                    const int y = ord == 0 ? qb : qa;
                    const int fx = I.df[static_cast<std::size_t>(x)];
                    const int fy = I.df[static_cast<std::size_t>(y)];
                    if (!((fx == f0 && fy == f1) || (fx == f1 && fy == f0))) {
                        continue;
                    }
                    const std::int8_t ox_end = (fx == f0) ? 0 : 1;
                    const std::int8_t oy_end = static_cast<std::int8_t>(1 - ox_end);
                    const PatOp P0{kPick, 0}, E0{kEmpty, 0}, P1{kPick, 1}, E1{kEmpty, 1}, D0{kDeliv, 0}, D1{kDeliv, 1};
                    if (all_skeletons) {
                        add(1, {{P0, E0, P1, E1, D0, D1}}, x, y, ox_end, oy_end);
                    }
                    add(3, {{P0, P1, E0, E1, D0, D1}}, x, y, ox_end, oy_end);
                    if (I.pb[static_cast<std::size_t>(p0)] != I.pb[static_cast<std::size_t>(p1)]) {
                        add(5, {{P0, P1, E1, E0, D0, D1}}, x, y, ox_end, oy_end);
                    }
                    if (all_skeletons) {
                        if (fx == f0 && fy == f1) {
                            add(2, {{P0, E0, P1, D0, E1, D1}}, x, y, 0, 1);
                            add(4, {{P0, P1, E0, D0, E1, D1}}, x, y, 0, 1);
                        }
                        if (fx == f1 && fy == f0) {
                            add(6, {{P0, P1, E1, D0, E0, D1}}, x, y, 1, 0);
                        }
                    }
                }
                struct Conc { int h, g; double d; };
                std::vector<Conc> concs;
                for (int ord = 0; ord < 2; ++ord) {
                    const int x = ord == 0 ? qa : qb;
                    const int y = ord == 0 ? qb : qa;
                    if (I.df[static_cast<std::size_t>(x)] != f0 || I.df[static_cast<std::size_t>(y)] != f1) {
                        continue;
                    }
                    const double d0 = I.single_d(p0, x);
                    const double d1 = I.single_d(p1, y);
                    concs.push_back(Conc{I.pa[static_cast<std::size_t>(p0)], I.dw[static_cast<std::size_t>(y)],
                        d0 + I.tt(I.dw[static_cast<std::size_t>(x)], I.pa[static_cast<std::size_t>(p1)]) + d1});
                    concs.push_back(Conc{I.pa[static_cast<std::size_t>(p1)], I.dw[static_cast<std::size_t>(x)],
                        d1 + I.tt(I.dw[static_cast<std::size_t>(y)], I.pa[static_cast<std::size_t>(p0)]) + d0});
                    if (qa == qb) {
                        break;
                    }
                }
                for (Composite& cp : cands) {
                    const int h = I.pa[static_cast<std::size_t>(p0)];
                    const int g = I.dw[static_cast<std::size_t>(cp.q1)];
                    bool dominated = false;
                    double best = kInf;
                    for (const Conc& s : concs) {
                        best = std::min(best, s.d);
                        if (s.h == h && s.g == g && s.d <= cp.d + kEps) {
                            dominated = true;
                        }
                    }
                    if (dominated) {
                        ++pool.dominated;
                        continue;
                    }
                    cp.delta = std::isfinite(best) ? best - cp.d : 0.0;
                    pool.comps.push_back(cp);
                    pool.comp_by_pair[pair_key(p0, p1, P)].push_back(
                        static_cast<int>(pool.comps.size()) - 1);
                }
            }
        }
    }
}

// ---------------------------------------------------------------- trips

struct Trip {
    bool composite = false;
    int comp = -1;
    int p0 = -1, p1 = -1, q0 = -1, q1 = -1;  // classes used
    double d = 0.0, c = 0.0, delta = 0.0;
    int ntraj = 1;
};

Trip make_single_trip(const Instance& I, int p, int q) {
    Trip t;
    t.composite = false;
    t.p0 = p;
    t.q0 = q;
    t.d = I.single_d(p, q);
    t.c = I.single_c(p, q);
    t.ntraj = 1;
    return t;
}

Trip make_composite_trip(const PatternPool& pool, int idx) {
    const Composite& cp = pool.comps[static_cast<std::size_t>(idx)];
    Trip t;
    t.composite = true;
    t.comp = idx;
    t.p0 = cp.p0; t.p1 = cp.p1; t.q0 = cp.q0; t.q1 = cp.q1;
    t.d = cp.d; t.c = cp.c; t.delta = cp.delta;
    t.ntraj = 2;
    return t;
}

std::vector<Op> trip_ops(const PatternPool& pool, const Trip& tr, int& next_traj) {
    std::vector<Op> ops;
    if (!tr.composite) {
        const int id = next_traj++;
        ops.push_back(Op{kPick, tr.p0, id});
        ops.push_back(Op{kEmpty, tr.p0, id});
        ops.push_back(Op{kDeliv, tr.q0, id});
        return ops;
    }
    const Composite& cp = pool.comps[static_cast<std::size_t>(tr.comp)];
    const int ids[2] = {next_traj, next_traj + 1};
    next_traj += 2;
    for (const PatOp& po : cp.ops) {
        if (po.kind == kDeliv) {
            ops.push_back(Op{kDeliv, po.slot ? cp.q1 : cp.q0, ids[cp.downer[static_cast<std::size_t>(po.slot)]]});
        } else {
            ops.push_back(Op{po.kind, po.slot ? cp.p1 : cp.p0, ids[po.slot]});
        }
    }
    return ops;
}

struct Counts {
    std::vector<int> u;  // remaining pickups per pickup class
    std::vector<int> z;  // remaining deliveries per delivery class
};

bool trip_valid(const Trip& t, const Counts& k) {
    if (!t.composite) {
        return k.u[static_cast<std::size_t>(t.p0)] >= 1 && k.z[static_cast<std::size_t>(t.q0)] >= 1;
    }
    if (t.p0 == t.p1) {
        if (k.u[static_cast<std::size_t>(t.p0)] < 2) return false;
    } else if (k.u[static_cast<std::size_t>(t.p0)] < 1 || k.u[static_cast<std::size_t>(t.p1)] < 1) {
        return false;
    }
    if (t.q0 == t.q1) {
        return k.z[static_cast<std::size_t>(t.q0)] >= 2;
    }
    return k.z[static_cast<std::size_t>(t.q0)] >= 1 && k.z[static_cast<std::size_t>(t.q1)] >= 1;
}

void consume_trip(const Trip& t, Counts& k) {
    --k.u[static_cast<std::size_t>(t.p0)];
    --k.z[static_cast<std::size_t>(t.q0)];
    if (t.composite) {
        --k.u[static_cast<std::size_t>(t.p1)];
        --k.z[static_cast<std::size_t>(t.q1)];
    }
}

// ---------------------------------------------------------------- assignment

// Square min-cost assignment (Hungarian, O(n^3)). Returns col for each row.
std::vector<int> hungarian(const std::vector<std::vector<double>>& a) {
    const int n = static_cast<int>(a.size());
    if (n == 0) return {};
    std::vector<double> u(static_cast<std::size_t>(n + 1), 0.0), v(static_cast<std::size_t>(n + 1), 0.0);
    std::vector<int> p(static_cast<std::size_t>(n + 1), 0), way(static_cast<std::size_t>(n + 1), 0);
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(static_cast<std::size_t>(n + 1), kInf);
        std::vector<char> used(static_cast<std::size_t>(n + 1), 0);
        do {
            used[static_cast<std::size_t>(j0)] = 1;
            const int i0 = p[static_cast<std::size_t>(j0)];
            double delta = kInf;
            int j1 = 0;
            for (int j = 1; j <= n; ++j) {
                if (used[static_cast<std::size_t>(j)]) continue;
                const double cur = a[static_cast<std::size_t>(i0 - 1)][static_cast<std::size_t>(j - 1)] -
                    u[static_cast<std::size_t>(i0)] - v[static_cast<std::size_t>(j)];
                if (cur < minv[static_cast<std::size_t>(j)]) {
                    minv[static_cast<std::size_t>(j)] = cur;
                    way[static_cast<std::size_t>(j)] = j0;
                }
                if (minv[static_cast<std::size_t>(j)] < delta) {
                    delta = minv[static_cast<std::size_t>(j)];
                    j1 = j;
                }
            }
            for (int j = 0; j <= n; ++j) {
                if (used[static_cast<std::size_t>(j)]) {
                    u[static_cast<std::size_t>(p[static_cast<std::size_t>(j)])] += delta;
                    v[static_cast<std::size_t>(j)] -= delta;
                } else {
                    minv[static_cast<std::size_t>(j)] -= delta;
                }
            }
            j0 = j1;
        } while (p[static_cast<std::size_t>(j0)] != 0);
        do {
            const int j1 = way[static_cast<std::size_t>(j0)];
            p[static_cast<std::size_t>(j0)] = p[static_cast<std::size_t>(j1)];
            j0 = j1;
        } while (j0 != 0);
    }
    std::vector<int> ans(static_cast<std::size_t>(n), -1);
    for (int j = 1; j <= n; ++j) {
        if (p[static_cast<std::size_t>(j)] > 0) {
            ans[static_cast<std::size_t>(p[static_cast<std::size_t>(j)] - 1)] = j - 1;
        }
    }
    return ans;
}

// mode: 0 greedy (min-cost assignment), 1 random matching, 2 sequential roulette.
std::vector<Trip> assign_singles(const Instance& I, const Counts& k, int mode, std::mt19937_64& rng,
                                 bool by_cost = false) {
    std::vector<Trip> out;
    std::vector<std::vector<int>> ulist(static_cast<std::size_t>(I.ntype)), zlist(static_cast<std::size_t>(I.ntype));
    for (std::size_t p = 0; p < k.u.size(); ++p) {
        for (int m = 0; m < k.u[p]; ++m) ulist[static_cast<std::size_t>(I.pf[p])].push_back(static_cast<int>(p));
    }
    for (std::size_t q = 0; q < k.z.size(); ++q) {
        for (int m = 0; m < k.z[q]; ++m) zlist[static_cast<std::size_t>(I.df[q])].push_back(static_cast<int>(q));
    }
    for (int f = 0; f < I.ntype; ++f) {
        std::vector<int>& U = ulist[static_cast<std::size_t>(f)];
        std::vector<int>& Z = zlist[static_cast<std::size_t>(f)];
        if (U.size() != Z.size()) {
            throw std::runtime_error("PALNS assignment lost the type balance.");
        }
        if (U.empty()) continue;
        if (mode == 0) {
            std::vector<std::vector<double>> cost(U.size(), std::vector<double>(Z.size(), 0.0));
            for (std::size_t i = 0; i < U.size(); ++i) {
                for (std::size_t j = 0; j < Z.size(); ++j) {
                    cost[i][j] = by_cost ? I.single_c(U[i], Z[j]) + 1e-6 * I.single_d(U[i], Z[j])
                                         : I.single_d(U[i], Z[j]) + 1e-6 * I.single_c(U[i], Z[j]);
                }
            }
            const std::vector<int> col = hungarian(cost);
            for (std::size_t i = 0; i < U.size(); ++i) {
                out.push_back(make_single_trip(I, U[i], Z[static_cast<std::size_t>(col[i])]));
            }
        } else {
            std::shuffle(U.begin(), U.end(), rng);
            std::vector<int> rem = Z;
            for (int p : U) {
                std::size_t pick = 0;
                if (mode == 1) {
                    std::uniform_int_distribution<std::size_t> dist(0, rem.size() - 1U);
                    pick = dist(rng);
                } else {
                    std::vector<double> w(rem.size(), 0.0);
                    double total = 0.0;
                    for (std::size_t j = 0; j < rem.size(); ++j) {
                        w[j] = by_cost ? 1.0 / (1.0 + I.single_c(p, rem[j]))
                                       : 1.0 / std::max(1e-9, I.single_d(p, rem[j]));
                        total += w[j];
                    }
                    std::uniform_real_distribution<double> dist(0.0, total);
                    double r = dist(rng);
                    pick = rem.size() - 1U;
                    for (std::size_t j = 0; j < rem.size(); ++j) {
                        r -= w[j];
                        if (r <= 0.0) { pick = j; break; }
                    }
                }
                out.push_back(make_single_trip(I, p, rem[pick]));
                rem.erase(rem.begin() + static_cast<std::ptrdiff_t>(pick));
            }
        }
    }
    return out;
}

// ---------------------------------------------------------------- insertion

struct InsEval {
    int route = -1;
    bool split = false;
    int i = -1, j = -1, k = -1;  // contiguous: i is the gap; split: P,E,D gaps
    double dV = kInf, dd = kInf, dc = kInf;
};

bool ins_better(const InsEval& a, const InsEval& b, bool phase2, double lambda) {
    if (a.route < 0) return false;
    if (b.route < 0) return true;
    if (phase2) {
        const double ka = a.dc + lambda * a.dV;
        const double kb = b.dc + lambda * b.dV;
        if (std::fabs(ka - kb) > kEps) return ka < kb;
        return a.dd < b.dd - kEps;
    }
    if (std::fabs(a.dV - b.dV) > kEps) return a.dV < b.dV;
    if (std::fabs(a.dd - b.dd) > kEps) return a.dd < b.dd;
    return a.dc < b.dc - kEps;
}

struct BlockInfo {
    std::vector<Op> ops;
    double d = 0.0, c = 0.0;
    int h = -1, g = -1;
    int lam = 0;
    bool single = false;
    int a = -1, b = -1, w = -1;  // single locations
};

BlockInfo make_block(const Instance& I, std::vector<Op> ops) {
    BlockInfo B;
    B.ops = std::move(ops);
    eval_ops(I, B.ops, B.d, B.c);
    B.h = I.loc(B.ops.front());
    B.g = I.loc(B.ops.back());
    int load = 0;
    for (const Op& o : B.ops) {
        load += load_delta(o);
        B.lam = std::max(B.lam, load);
    }
    B.single = B.ops.size() == 3U;
    if (B.single) {
        B.a = I.loc(B.ops[0]);
        B.b = I.loc(B.ops[1]);
        B.w = I.loc(B.ops[2]);
    }
    return B;
}

struct RouteCtx {
    std::vector<int> L;
    std::vector<int> lo;
};

void make_ctx(const Instance& I, const Route& r, RouteCtx& ctx) {
    route_loads(r.ops, ctx.L);
    ctx.lo.resize(r.ops.size());
    for (std::size_t i = 0; i < r.ops.size(); ++i) {
        ctx.lo[i] = I.loc(r.ops[i]);
    }
}

// Cost of inserting a contiguous location chain (first x, last y, internal
// duration/cost) at gap g.
inline void chain_at(const Instance& I, const RouteCtx& ctx, int g, int x, int y,
                     double inner_d, double inner_c, double& dd, double& dc) {
    const int m = static_cast<int>(ctx.lo.size());
    const int prev = g > 0 ? ctx.lo[static_cast<std::size_t>(g - 1)] : -1;
    const int next = g < m ? ctx.lo[static_cast<std::size_t>(g)] : -1;
    dd = inner_d + I.lt(prev, x) + I.lt(y, next) - I.lt(prev, next);
    dc = inner_c + I.lc(prev, x) + I.lc(y, next) - I.lc(prev, next);
}

void exact_split(const Instance& I, const RouteCtx& ctx, const BlockInfo& B,
                 int i, int j, int k, double& dd, double& dc) {
    double d1 = 0, c1 = 0, d2 = 0, c2 = 0, d3 = 0, c3 = 0;
    if (i == j && j == k) {
        chain_at(I, ctx, i, B.a, B.w, B.d, B.c, dd, dc);
        return;
    }
    if (i == j) {
        chain_at(I, ctx, i, B.a, B.b, I.sp + I.se + I.tt(B.a, B.b), I.cc(B.a, B.b), d1, c1);
        chain_at(I, ctx, k, B.w, B.w, I.sd, 0.0, d2, c2);
        dd = d1 + d2; dc = c1 + c2;
        return;
    }
    if (j == k) {
        chain_at(I, ctx, i, B.a, B.a, I.sp, 0.0, d1, c1);
        chain_at(I, ctx, j, B.b, B.w, I.se + I.sd + I.tt(B.b, B.w), I.cc(B.b, B.w), d2, c2);
        dd = d1 + d2; dc = c1 + c2;
        return;
    }
    chain_at(I, ctx, i, B.a, B.a, I.sp, 0.0, d1, c1);
    chain_at(I, ctx, j, B.b, B.b, I.se, 0.0, d2, c2);
    chain_at(I, ctx, k, B.w, B.w, I.sd, 0.0, d3, c3);
    dd = d1 + d2 + d3; dc = c1 + c2 + c3;
}

// Best split position (P at i, E at j, D at k) by duration (use_cost=false)
// or by cost (use_cost=true) inside windows with L <= 1.
void best_split_by(const Instance& I, const RouteCtx& ctx, const BlockInfo& B, bool use_cost,
                   int& bi, int& bj, int& bk, double& bv) {
    const int m = static_cast<int>(ctx.lo.size());
    bi = bj = bk = -1;
    bv = kInf;
    auto val = [&](double d, double c) { return use_cost ? c : d; };
    int g = 0;
    while (g <= m) {
        if (ctx.L[static_cast<std::size_t>(g)] > 1) { ++g; continue; }
        const int lo_g = g;
        while (g <= m && ctx.L[static_cast<std::size_t>(g)] <= 1) ++g;
        const int hi_g = g - 1;
        const int len = hi_g - lo_g + 1;
        std::vector<double> dP(static_cast<std::size_t>(len)), dE(static_cast<std::size_t>(len)), dD(static_cast<std::size_t>(len));
        std::vector<double> pe(static_cast<std::size_t>(len)), ed(static_cast<std::size_t>(len)), ped(static_cast<std::size_t>(len));
        for (int x = 0; x < len; ++x) {
            const int gg = lo_g + x;
            double d, c;
            chain_at(I, ctx, gg, B.a, B.a, I.sp, 0.0, d, c); dP[static_cast<std::size_t>(x)] = val(d, c);
            chain_at(I, ctx, gg, B.b, B.b, I.se, 0.0, d, c); dE[static_cast<std::size_t>(x)] = val(d, c);
            chain_at(I, ctx, gg, B.w, B.w, I.sd, 0.0, d, c); dD[static_cast<std::size_t>(x)] = val(d, c);
            chain_at(I, ctx, gg, B.a, B.b, I.sp + I.se + I.tt(B.a, B.b), I.cc(B.a, B.b), d, c); pe[static_cast<std::size_t>(x)] = val(d, c);
            chain_at(I, ctx, gg, B.b, B.w, I.se + I.sd + I.tt(B.b, B.w), I.cc(B.b, B.w), d, c); ed[static_cast<std::size_t>(x)] = val(d, c);
            chain_at(I, ctx, gg, B.a, B.w, B.d, B.c, d, c); ped[static_cast<std::size_t>(x)] = val(d, c);
        }
        std::vector<double> preV(static_cast<std::size_t>(len), kInf), sufV(static_cast<std::size_t>(len), kInf);
        std::vector<int> preI(static_cast<std::size_t>(len), -1), sufI(static_cast<std::size_t>(len), -1);
        for (int x = 0; x < len; ++x) {
            if (x > 0) { preV[static_cast<std::size_t>(x)] = preV[static_cast<std::size_t>(x - 1)]; preI[static_cast<std::size_t>(x)] = preI[static_cast<std::size_t>(x - 1)]; }
            if (dP[static_cast<std::size_t>(x)] < preV[static_cast<std::size_t>(x)]) { preV[static_cast<std::size_t>(x)] = dP[static_cast<std::size_t>(x)]; preI[static_cast<std::size_t>(x)] = x; }
        }
        for (int x = len - 1; x >= 0; --x) {
            if (x + 1 < len) { sufV[static_cast<std::size_t>(x)] = sufV[static_cast<std::size_t>(x + 1)]; sufI[static_cast<std::size_t>(x)] = sufI[static_cast<std::size_t>(x + 1)]; }
            if (dD[static_cast<std::size_t>(x)] < sufV[static_cast<std::size_t>(x)]) { sufV[static_cast<std::size_t>(x)] = dD[static_cast<std::size_t>(x)]; sufI[static_cast<std::size_t>(x)] = x; }
        }
        auto consider = [&](double v, int i, int j, int k) {
            if (v < bv - 1e-12) { bv = v; bi = lo_g + i; bj = lo_g + j; bk = lo_g + k; }
        };
        for (int x = 0; x < len; ++x) {
            consider(ped[static_cast<std::size_t>(x)], x, x, x);
            if (x + 1 < len) consider(pe[static_cast<std::size_t>(x)] + sufV[static_cast<std::size_t>(x + 1)], x, x, sufI[static_cast<std::size_t>(x + 1)]);
            if (x > 0) consider(preV[static_cast<std::size_t>(x - 1)] + ed[static_cast<std::size_t>(x)], preI[static_cast<std::size_t>(x - 1)], x, x);
            if (x > 0 && x + 1 < len) {
                consider(preV[static_cast<std::size_t>(x - 1)] + dE[static_cast<std::size_t>(x)] + sufV[static_cast<std::size_t>(x + 1)],
                         preI[static_cast<std::size_t>(x - 1)], x, sufI[static_cast<std::size_t>(x + 1)]);
            }
        }
    }
}

// Scratch arrays of the split-position search (version 2), reused between calls.
struct SplitScratch {
    std::vector<double> dP[2], dE[2], dD[2], pe[2], ed[2], ped[2];  // [0] duration, [1] cost
    std::vector<double> preV, sufV;
    std::vector<int> preI, sufI;
    void size(std::size_t n) {
        for (int k = 0; k < 2; ++k) {
            dP[k].resize(n); dE[k].resize(n); dD[k].resize(n);
            pe[k].resize(n); ed[k].resize(n); ped[k].resize(n);
        }
        preV.resize(n); sufV.resize(n); preI.resize(n); sufI.resize(n);
    }
};

// Best split positions by duration ([0]) and, when both is set, by cost ([1]):
// the result of best_split_by for use_cost = false / true, computed in one pass
// over the gaps (every chain is evaluated once for both criteria).
void best_split_both(const Instance& I, const RouteCtx& ctx, const BlockInfo& B, bool both,
                     SplitScratch& sc, int bi[2], int bj[2], int bk[2]) {
    const int m = static_cast<int>(ctx.lo.size());
    double bv[2] = {kInf, kInf};
    for (int k = 0; k < 2; ++k) bi[k] = bj[k] = bk[k] = -1;
    const int npass = both ? 2 : 1;
    int g = 0;
    while (g <= m) {
        if (ctx.L[static_cast<std::size_t>(g)] > 1) { ++g; continue; }
        const int lo_g = g;
        while (g <= m && ctx.L[static_cast<std::size_t>(g)] <= 1) ++g;
        const int len = g - lo_g;
        sc.size(static_cast<std::size_t>(len));
        for (int x = 0; x < len; ++x) {
            const int gg = lo_g + x;
            const std::size_t xs = static_cast<std::size_t>(x);
            double d, c;
            chain_at(I, ctx, gg, B.a, B.a, I.sp, 0.0, d, c); sc.dP[0][xs] = d; sc.dP[1][xs] = c;
            chain_at(I, ctx, gg, B.b, B.b, I.se, 0.0, d, c); sc.dE[0][xs] = d; sc.dE[1][xs] = c;
            chain_at(I, ctx, gg, B.w, B.w, I.sd, 0.0, d, c); sc.dD[0][xs] = d; sc.dD[1][xs] = c;
            chain_at(I, ctx, gg, B.a, B.b, I.sp + I.se + I.tt(B.a, B.b), I.cc(B.a, B.b), d, c); sc.pe[0][xs] = d; sc.pe[1][xs] = c;
            chain_at(I, ctx, gg, B.b, B.w, I.se + I.sd + I.tt(B.b, B.w), I.cc(B.b, B.w), d, c); sc.ed[0][xs] = d; sc.ed[1][xs] = c;
            chain_at(I, ctx, gg, B.a, B.w, B.d, B.c, d, c); sc.ped[0][xs] = d; sc.ped[1][xs] = c;
        }
        for (int k = 0; k < npass; ++k) {
            const std::vector<double>& dP = sc.dP[k];
            const std::vector<double>& dE = sc.dE[k];
            const std::vector<double>& dD = sc.dD[k];
            const std::vector<double>& pe = sc.pe[k];
            const std::vector<double>& ed = sc.ed[k];
            const std::vector<double>& ped = sc.ped[k];
            for (int x = 0; x < len; ++x) {
                const std::size_t xs = static_cast<std::size_t>(x);
                sc.preV[xs] = kInf; sc.preI[xs] = -1;
                if (x > 0) { sc.preV[xs] = sc.preV[xs - 1U]; sc.preI[xs] = sc.preI[xs - 1U]; }
                if (dP[xs] < sc.preV[xs]) { sc.preV[xs] = dP[xs]; sc.preI[xs] = x; }
            }
            for (int x = len - 1; x >= 0; --x) {
                const std::size_t xs = static_cast<std::size_t>(x);
                sc.sufV[xs] = kInf; sc.sufI[xs] = -1;
                if (x + 1 < len) { sc.sufV[xs] = sc.sufV[xs + 1U]; sc.sufI[xs] = sc.sufI[xs + 1U]; }
                if (dD[xs] < sc.sufV[xs]) { sc.sufV[xs] = dD[xs]; sc.sufI[xs] = x; }
            }
            auto consider = [&](double v, int i, int j, int kk) {
                if (v < bv[k] - 1e-12) { bv[k] = v; bi[k] = lo_g + i; bj[k] = lo_g + j; bk[k] = lo_g + kk; }
            };
            for (int x = 0; x < len; ++x) {
                const std::size_t xs = static_cast<std::size_t>(x);
                consider(ped[xs], x, x, x);
                if (x + 1 < len) consider(pe[xs] + sc.sufV[xs + 1U], x, x, sc.sufI[xs + 1U]);
                if (x > 0) consider(sc.preV[xs - 1U] + ed[xs], sc.preI[xs - 1U], x, x);
                if (x > 0 && x + 1 < len) {
                    consider(sc.preV[xs - 1U] + dE[xs] + sc.sufV[xs + 1U], sc.preI[xs - 1U], x, sc.sufI[xs + 1U]);
                }
            }
        }
    }
}

// Best insertion of block B into one route: the per-route step of
// eval_insertion, with the same candidates and the same tie-breaking. ctx is
// the context of the route (make_ctx).
InsEval eval_insertion_route(const Instance& I, const Route& route, const RouteCtx& ctx, int r,
                             const BlockInfo& B, bool allow_split, bool phase2, double lambda,
                             SplitScratch& sc) {
    const int m = static_cast<int>(route.ops.size());
    InsEval rbest;
    auto offer = [&](InsEval e) {
        e.route = r;
        e.dV = excess(route.d + e.dd, I.T) - excess(route.d, I.T);
        if (ins_better(e, rbest, phase2, lambda)) rbest = e;
    };
    for (int g = 0; g <= m; ++g) {
        if (ctx.L[static_cast<std::size_t>(g)] + B.lam > 2) continue;
        InsEval e;
        e.split = false;
        e.i = g;
        chain_at(I, ctx, g, B.h, B.g, B.d, B.c, e.dd, e.dc);
        offer(e);
    }
    if (allow_split && B.single) {
        int bi[2], bj[2], bk[2];
        best_split_both(I, ctx, B, phase2, sc, bi, bj, bk);
        for (int pass = 0; pass < (phase2 ? 2 : 1); ++pass) {
            if (bi[pass] < 0) continue;
            InsEval e;
            e.split = true;
            e.i = bi[pass]; e.j = bj[pass]; e.k = bk[pass];
            exact_split(I, ctx, B, bi[pass], bj[pass], bk[pass], e.dd, e.dc);
            offer(e);
        }
    }
    return rbest;
}

// Evaluates the best insertion of block B into every route. per_route (if
// given) receives the best position per route (for regret). Same result as in
// version 1; the per-route step is eval_insertion_route.
InsEval eval_insertion(const Instance& I, const Solution& s, const BlockInfo& B, bool allow_split,
                       bool phase2, double lambda, std::vector<InsEval>* per_route) {
    static thread_local RouteCtx ctx;
    static thread_local SplitScratch sc;
    InsEval best;
    if (per_route != nullptr) per_route->assign(s.routes.size(), InsEval{});
    for (std::size_t r = 0; r < s.routes.size(); ++r) {
        make_ctx(I, s.routes[r], ctx);
        const InsEval rbest = eval_insertion_route(I, s.routes[r], ctx, static_cast<int>(r), B, allow_split,
                                                   phase2, lambda, sc);
        if (per_route != nullptr) (*per_route)[r] = rbest;
        if (ins_better(rbest, best, phase2, lambda)) best = rbest;
    }
    return best;
}

void apply_insertion(const Instance& I, Solution& s, const BlockInfo& B, const InsEval& e) {
    Route& route = s.routes[static_cast<std::size_t>(e.route)];
    const double old_d = route.d;
    const double old_c = route.c;
    std::vector<Op> res;
    res.reserve(route.ops.size() + B.ops.size());
    const int m = static_cast<int>(route.ops.size());
    for (int g = 0; g <= m; ++g) {
        if (!e.split) {
            if (g == e.i) res.insert(res.end(), B.ops.begin(), B.ops.end());
        } else {
            if (g == e.i) res.push_back(B.ops[0]);
            if (g == e.j) res.push_back(B.ops[1]);
            if (g == e.k) res.push_back(B.ops[2]);
        }
        if (g < m) res.push_back(route.ops[static_cast<std::size_t>(g)]);
    }
    route.ops.swap(res);
    refresh_route(I, route);
    s.V += excess(route.d, I.T) - excess(old_d, I.T);
    s.D += route.d - old_d;
    s.C += route.c - old_c;
}

// ---------------------------------------------------------------- helpers

struct TrajLoc {
    int route = -1;
    int pi = -1, ei = -1, di = -1;
    int p = -1, q = -1;
};

std::unordered_map<int, TrajLoc> collect_trajs(const Solution& s) {
    std::unordered_map<int, TrajLoc> out;
    for (std::size_t r = 0; r < s.routes.size(); ++r) {
        const auto& ops = s.routes[r].ops;
        for (std::size_t i = 0; i < ops.size(); ++i) {
            TrajLoc& t = out[ops[i].traj];
            t.route = static_cast<int>(r);
            if (ops[i].kind == kPick) { t.pi = static_cast<int>(i); t.p = ops[i].cls; }
            else if (ops[i].kind == kEmpty) { t.ei = static_cast<int>(i); }
            else { t.di = static_cast<int>(i); t.q = ops[i].cls; }
        }
    }
    return out;
}

void remove_trajs(const Instance& I, Solution& s, const std::unordered_set<int>& rem, Counts& k) {
    for (Route& route : s.routes) {
        bool touched = false;
        for (const Op& o : route.ops) {
            if (rem.count(o.traj) != 0U) { touched = true; break; }
        }
        if (!touched) continue;
        std::vector<Op> kept;
        kept.reserve(route.ops.size());
        for (const Op& o : route.ops) {
            if (rem.count(o.traj) != 0U) {
                if (o.kind == kPick) ++k.u[static_cast<std::size_t>(o.cls)];
                else if (o.kind == kDeliv) ++k.z[static_cast<std::size_t>(o.cls)];
                continue;
            }
            kept.push_back(o);
        }
        route.ops.swap(kept);
        refresh_route(I, route);
    }
    refresh_totals(I, s);
}

void route_blocks(const std::vector<Op>& ops, std::vector<std::pair<int, int>>& blocks) {
    blocks.clear();
    int load = 0;
    int start = 0;
    for (std::size_t i = 0; i < ops.size(); ++i) {
        load += load_delta(ops[i]);
        if (load == 0) {
            blocks.push_back({start, static_cast<int>(i) + 1});
            start = static_cast<int>(i) + 1;
        }
    }
}

std::uint64_t splitmix(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

std::uint64_t solution_hash(const Solution& s) {
    std::uint64_t h = 0;
    for (const Route& r : s.routes) {
        std::uint64_t prev = 1;
        for (const Op& o : r.ops) {
            const std::uint64_t lab = (static_cast<std::uint64_t>(o.kind + 1) << 32) |
                static_cast<std::uint32_t>(o.cls);
            h += splitmix(prev * 1000003ULL + lab);
            prev = lab;
        }
        h += splitmix(prev * 1000003ULL + 2ULL);
    }
    return h;
}

template <typename T>
int roulette_index(const std::vector<T>& weights, std::mt19937_64& rng) {
    double total = 0.0;
    for (const T& w : weights) total += static_cast<double>(w);
    if (total <= 0.0) {
        std::uniform_int_distribution<int> dist(0, static_cast<int>(weights.size()) - 1);
        return dist(rng);
    }
    std::uniform_real_distribution<double> dist(0.0, total);
    double r = dist(rng);
    for (std::size_t i = 0; i < weights.size(); ++i) {
        r -= static_cast<double>(weights[i]);
        if (r <= 0.0) return static_cast<int>(i);
    }
    return static_cast<int>(weights.size()) - 1;
}

// Version 2 switches by name (see PalnsV2Options). Returns false for an unknown key;
// std::sto* throw on a malformed value.
bool set_v2_option(PalnsV2Options& o, const std::string& key, const std::string& value) {
    if (key == "sel_cost") o.sel_cost = std::stoi(value) != 0;
    else if (key == "fast_relocate") o.fast_relocate = std::stoi(value) != 0;
    else if (key == "profile") o.profile = std::stoi(value) != 0;
    else if (key == "diag") o.diag = std::stoi(value) != 0;
    else if (key == "temp_spread") o.temp_spread = std::stod(value);
    else if (key == "sel_cost_threads") o.sel_cost_threads = std::stoi(value);
    else if (key == "focus_period") o.focus_period = std::stoi(value);
    else if (key == "focus_routes") o.focus_routes = std::stoi(value);
    else if (key == "focus_len") o.focus_len = std::stoi(value);
    else if (key == "focus_kmax_ratio") o.focus_kmax_ratio = std::stod(value);
    else if (key == "focus_rebuild") o.focus_rebuild = std::stod(value);
    else if (key == "focus_deep") o.focus_deep = std::stoi(value);
    else if (key == "focus_start_ratio") o.focus_start_ratio = std::stod(value);
    else if (key == "lam_lo") o.lam_lo = std::stod(value);
    else if (key == "lam_hi") o.lam_hi = std::stod(value);
    else if (key == "lam_step") o.lam_step = std::stod(value);
    else if (key == "lam_dstep") o.lam_dstep = std::stod(value);
    else if (key == "lam_win") o.lam_win = std::stoi(value);
    else if (key == "lam_sig") o.lam_sig = std::stoi(value);
    else if (key == "lam_back") o.lam_back = std::stoi(value);
    else if (key == "lam_to") o.lam_to = std::stoi(value);
    else if (key == "lam_tlo") o.lam_tlo = std::stod(value);
    else if (key == "lam_thi") o.lam_thi = std::stod(value);
    else if (key == "p1_join") o.p1_join = std::stoi(value);
    else if (key == "indep_adopt") o.indep_adopt = std::stoi(value);
    else if (key == "reheat_max") o.reheat_max = std::stod(value);
    else if (key == "reheat_step") o.reheat_step = std::stod(value);
    else if (key == "reheat_after") o.reheat_after = std::stoi(value);
    else if (key == "rescue_p2") o.rescue_p2 = std::stoi(value) != 0;
    else if (key == "pool_ratio") o.pool_ratio = std::stod(value);
    else if (key == "pool_file") o.pool_file = value;
    else if (key == "pool_delta") o.pool_delta = std::stod(value);
    else if (key == "recomb_tl") o.recomb_tl = std::stod(value);
    else if (key == "recomb_delta") o.recomb_delta = std::stod(value);
    else if (key == "recomb_threads") o.recomb_threads = std::stoi(value);
    else return false;
    return true;
}

}  // namespace

// ---------------------------------------------------------------- solver

namespace {
// Epoch barrier of the portfolio. Threads run sync_iters iterations, then
// wait here; the last arriving thread runs on_barrier (merge of the thread
// bests in thread order, stop rule) while all others are blocked. Between two
// barriers the shared best snapshot below is read-only, so the search is
// deterministic for a fixed seed.
struct SyncHub {
    std::mutex m;
    std::condition_variable cv;
    int participants = 0;
    int arrived = 0;
    long long generation = 0;
    bool stop = false;
    bool rescue = false;    // Phase I limit reached: threads switch to rescue
    bool rhas = false;      // rescue: best infeasible solution over the threads, (V, D)
    Solution rsol;
    std::function<void()> on_barrier;
    bool has = false;       // shared best snapshot
    long long first_i = -1; // iteration of the barrier at which it first existed
    double cost = 0.0;
    Solution solution;
};

bool hub_arrive(SyncHub& h, double& wait_seconds) {
    const auto w0 = Clock::now();
    std::unique_lock<std::mutex> lk(h.m);
    const long long gen = h.generation;
    if (++h.arrived >= h.participants) {
        if (h.on_barrier) h.on_barrier();
        h.arrived = 0;
        ++h.generation;
        h.cv.notify_all();
    } else {
        h.cv.wait(lk, [&]() { return h.generation != gen; });
    }
    wait_seconds += std::chrono::duration<double>(Clock::now() - w0).count();
    return h.stop;
}

void hub_leave(SyncHub& h) {
    std::unique_lock<std::mutex> lk(h.m);
    --h.participants;
    if (h.participants > 0 && h.arrived > 0 && h.arrived >= h.participants) {
        if (h.on_barrier) h.on_barrier();
        h.arrived = 0;
        ++h.generation;
        h.cv.notify_all();
    }
}
}  // namespace

struct PalnsV2IncumbentSolver::Impl {
    const SPDPData& data;
    PalnsOptions opt;
    PalnsV2Options v2;  // version 2 switches (all off: identical to version 1)
    Instance I;
    PatternPool pool;
    double pattern_seconds = 0.0;
    int next_traj = 0;
    int phase = 1;
    double lambda = 1.0;  // Phase II penalty on duration excess
    double eta_cur = 0.01;  // Phase I weight of the total duration
    Clock::time_point t0;
    Clock::time_point deadline;
    bool has_deadline = false;
    SyncHub* hub = nullptr;
    bool adopt_shared = true;  // odd portfolio threads keep searching independently
    long long gi = 0;          // iterations of this thread since the start (both phases)
    int rescue_switches = 0;   // 1 when this thread entered the rescue
    double rescue_start_V = -1.0;  // violation of the rescue start solution
    std::vector<std::size_t> rescue_nb;  // current rescue neighbourhood
    int rescue_nb_left = 0;              // iterations left on it
    bool rbest_has = false;              // rescue best of this thread, (V, D)
    Solution rbest;
    double wait_seconds = 0.0; // time spent at epoch barriers
    // Best feasible solution of this thread (merged at the next barrier).
    bool cand_has = false;
    double cand_cost = 0.0;
    long long cand_gi = 0;
    Solution cand;
    // Version 2: state of the focus blocks and counters (per thread).
    struct Episode {
        bool active = false;
        std::vector<std::size_t> nb;   // routes of the thread solution that the episode works on
        Solution cur, best;            // current and best feasible sub-solution (routes in the order of nb)
        double fcur = 0.0, theta0 = 0.0;
        long long len = 0, done = 0;
    };
    Episode ep;
    std::vector<long long> focus_last;  // episode number at which each route was last the seed
    long long focus_count = 0;
    long long stat_focus_iters = 0, stat_focus_episodes = 0, stat_focus_best = 0, stat_global_best = 0;

    void renumber(Solution& s) {
        std::unordered_map<int, int> remap;
        for (Route& r : s.routes) {
            for (Op& o : r.ops) {
                auto it = remap.find(o.traj);
                if (it == remap.end()) it = remap.emplace(o.traj, next_traj++).first;
                o.traj = it->second;
            }
        }
    }

    // Records a feasible solution of this thread; it reaches the shared best
    // only at the next epoch barrier.
    void publish(const Solution& s) {
        if (s.V > kEps) return;
        if (!cand_has || s.C < cand_cost - kEps) {
            cand_has = true;
            cand_cost = s.C;
            cand_gi = gi;
            cand = s;
        }
    }

    // Reads the shared best snapshot taken at the last barrier.
    bool adopt(Solution& out, double own_cost, bool force = false) {
        if (hub == nullptr || (!adopt_shared && !force)) return false;
        if (!hub->has || hub->cost >= own_cost - kEps) return false;
        out = hub->solution;
        renumber(out);
        return true;
    }

    Impl(const SPDPData& d, const PalnsOptions& o) : data(d), opt(o) {
        const auto start = Clock::now();
        build_instance();
        build_pattern_pool(I, opt, pool);
        pattern_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    }

    void build_instance() {
        I.n = static_cast<int>(data.requests.size());
        I.sp = data.time_pickup; I.se = data.time_empty; I.sd = data.time_delivery;
        I.T = data.time_limit; I.F = data.fixed_vehicle_cost;
        I.t = &data.time; I.c = &data.distance;
        std::map<int, int> type_id;
        for (const Request& r : data.requests) {
            if (type_id.find(r.container_type) == type_id.end()) {
                const int id = static_cast<int>(type_id.size());
                type_id[r.container_type] = id;
            }
        }
        I.ntype = static_cast<int>(type_id.size());
        std::map<std::array<int, 3>, int> pcls;
        std::map<std::array<int, 2>, int> dcls;
        I.req_pc.assign(static_cast<std::size_t>(I.n), -1);
        I.req_dc.assign(static_cast<std::size_t>(I.n), -1);
        for (int r = 0; r < I.n; ++r) {
            const Request& req = data.requests[static_cast<std::size_t>(r)];
            const int f = type_id[req.container_type];
            const std::array<int, 3> pk{{req.from_id, req.to_id, f}};
            auto it = pcls.find(pk);
            if (it == pcls.end()) {
                const int id = static_cast<int>(pcls.size());
                it = pcls.emplace(pk, id).first;
                I.pa.push_back(req.from_id); I.pb.push_back(req.to_id); I.pf.push_back(f);
                I.pm.push_back(0); I.pmem.emplace_back();
            }
            I.req_pc[static_cast<std::size_t>(r)] = it->second;
            ++I.pm[static_cast<std::size_t>(it->second)];
            I.pmem[static_cast<std::size_t>(it->second)].push_back(r);
            const std::array<int, 2> dk{{req.from_id, f}};
            auto jt = dcls.find(dk);
            if (jt == dcls.end()) {
                const int id = static_cast<int>(dcls.size());
                jt = dcls.emplace(dk, id).first;
                I.dw.push_back(req.from_id); I.df.push_back(f);
                I.dm.push_back(0); I.dmem.emplace_back();
            }
            I.req_dc[static_cast<std::size_t>(r)] = jt->second;
            ++I.dm[static_cast<std::size_t>(jt->second)];
            I.dmem[static_cast<std::size_t>(jt->second)].push_back(r);
        }
    }

    int npc() const { return static_cast<int>(I.pa.size()); }
    int ndc() const { return static_cast<int>(I.dw.size()); }

    double elapsed() const {
        return std::chrono::duration<double>(Clock::now() - t0).count();
    }
    bool timed_out() const {
        if (hub != nullptr && hub->stop) return true;
        return has_deadline && Clock::now() >= deadline;
    }

    double objective(const Solution& s) const {
        return phase == 1 ? s.V + eta_cur * s.D : s.C + lambda * s.V;
    }
    double original_cost(const Solution& s) const {
        return s.C + I.F * nonempty_route_count(s);
    }

    int select_comp(const std::vector<int>& list, int mode, std::mt19937_64& rng,
                    const std::vector<int>* zc) const {
        std::vector<int> ok;
        ok.reserve(list.size());
        for (int idx : list) {
            const Composite& cp = pool.comps[static_cast<std::size_t>(idx)];
            if (zc != nullptr) {
                if (cp.q0 == cp.q1) {
                    if ((*zc)[static_cast<std::size_t>(cp.q0)] < 2) continue;
                } else if ((*zc)[static_cast<std::size_t>(cp.q0)] < 1 || (*zc)[static_cast<std::size_t>(cp.q1)] < 1) {
                    continue;
                }
            }
            ok.push_back(idx);
        }
        if (ok.empty()) return -1;
        const bool by_cost = phase == 2 && v2.sel_cost;
        if (mode == 0) {
            int best = ok.front();
            for (int idx : ok) {
                const Composite& a = pool.comps[static_cast<std::size_t>(idx)];
                const Composite& b = pool.comps[static_cast<std::size_t>(best)];
                const double ka = by_cost ? a.c : a.d, kb = by_cost ? b.c : b.d;
                const double ta = by_cost ? a.d : a.c, tb = by_cost ? b.d : b.c;
                if (ka < kb - kEps || (std::fabs(ka - kb) <= kEps && ta < tb - kEps)) best = idx;
            }
            return best;
        }
        if (mode == 1) {
            std::uniform_int_distribution<std::size_t> dist(0, ok.size() - 1U);
            return ok[dist(rng)];
        }
        std::vector<double> w(ok.size());
        for (std::size_t i = 0; i < ok.size(); ++i) {
            const Composite& cpi = pool.comps[static_cast<std::size_t>(ok[i])];
            w[i] = by_cost ? 1.0 / (1.0 + cpi.c) : 1.0 / std::max(1e-9, cpi.d);
        }
        return ok[static_cast<std::size_t>(roulette_index(w, rng))];
    }

    // ------------------------------------------------------------ insertion helpers
    void insert_trip(Solution& s, const Trip& t, bool allow_split, bool phase2) {
        BlockInfo B = make_block(I, trip_ops(pool, t, next_traj));
        const InsEval e = eval_insertion(I, s, B, allow_split, phase2, lambda, nullptr);
        if (e.route < 0) {
            throw std::runtime_error("PALNS found no insertion position.");
        }
        if (v2.diag && phase == 2) {
            if (t.composite) ++diag.comp_ins; else ++diag.single_ins;
            if (e.split) ++diag.single_split;
        }
        apply_insertion(I, s, B, e);
    }

    // Regret/greedy ordered insertion of a candidate list. Returns nothing;
    // counts are consumed when non-null.
    void insert_dynamic(Solution& s, std::vector<Trip>& cand, Counts* cnt, bool regret,
                        bool allow_split, bool phase2) {
        std::vector<BlockInfo> blocks;
        blocks.reserve(cand.size());
        for (const Trip& t : cand) blocks.push_back(make_block(I, trip_ops(pool, t, next_traj)));
        std::vector<char> done(cand.size(), 0);
        // Version 2: the best position of a candidate in a route is recomputed
        // only for the route changed by the last insertion (same values as a
        // full re-evaluation, since the other routes did not change).
        const std::size_t R = s.routes.size();
        std::vector<std::vector<InsEval>> evals(cand.size());
        std::vector<RouteCtx> ctxs(R);
        for (std::size_t r = 0; r < R; ++r) make_ctx(I, s.routes[r], ctxs[r]);
        int dirty = -1;  // route changed by the last insertion
        while (true) {
            int pick = -1;
            InsEval pick_eval;
            double pk1 = kInf, pk2 = kInf, pk3 = kInf;  // ordering keys
            for (std::size_t x = 0; x < cand.size(); ++x) {
                if (done[x]) continue;
                if (cnt != nullptr && !trip_valid(cand[x], *cnt)) { done[x] = 1; continue; }
                std::vector<InsEval>& per_route = evals[x];
                if (per_route.empty()) {
                    per_route.resize(R);
                    for (std::size_t r = 0; r < R; ++r) {
                        per_route[r] = eval_insertion_route(I, s.routes[r], ctxs[r], static_cast<int>(r), blocks[x],
                                                            allow_split, phase2, lambda, ins_scratch);
                    }
                } else if (dirty >= 0) {
                    const std::size_t r = static_cast<std::size_t>(dirty);
                    per_route[r] = eval_insertion_route(I, s.routes[r], ctxs[r], dirty, blocks[x],
                                                        allow_split, phase2, lambda, ins_scratch);
                }
                InsEval e;
                for (std::size_t r = 0; r < R; ++r) {
                    if (ins_better(per_route[r], e, phase2, lambda)) e = per_route[r];
                }
                if (e.route < 0) continue;
                double k1, k2, k3;
                const double nt = static_cast<double>(cand[x].ntraj);
                if (!regret) {
                    k1 = phase2 ? (e.dc + lambda * e.dV) / nt : e.dV / nt;
                    k2 = e.dd / nt;
                    k3 = 0.0;
                } else {
                    InsEval second;
                    for (const InsEval& r : per_route) {
                        if (r.route < 0 || r.route == e.route) continue;
                        if (ins_better(r, second, phase2, lambda)) second = r;
                    }
                    double rv = kInf, rd = kInf;
                    if (second.route >= 0) {
                        rv = phase2 ? (second.dc + lambda * second.dV) - (e.dc + lambda * e.dV)
                                    : second.dV - e.dV;
                        rd = second.dd - e.dd;
                    }
                    k1 = -rv; k2 = -rd; k3 = e.dd;
                }
                bool better = false;
                if (pick < 0) better = true;
                else if (std::fabs(k1 - pk1) > kEps) better = k1 < pk1;
                else if (std::fabs(k2 - pk2) > kEps) better = k2 < pk2;
                else better = k3 < pk3 - kEps;
                if (better) { pick = static_cast<int>(x); pick_eval = e; pk1 = k1; pk2 = k2; pk3 = k3; }
            }
            if (pick < 0) break;
            if (v2.diag && phase == 2) {
                if (cand[static_cast<std::size_t>(pick)].composite) ++diag.comp_ins; else ++diag.single_ins;
                if (pick_eval.split) ++diag.single_split;
            }
            apply_insertion(I, s, blocks[static_cast<std::size_t>(pick)], pick_eval);
            if (cnt != nullptr) consume_trip(cand[static_cast<std::size_t>(pick)], *cnt);
            done[static_cast<std::size_t>(pick)] = 1;
            dirty = pick_eval.route;
            make_ctx(I, s.routes[static_cast<std::size_t>(dirty)], ctxs[static_cast<std::size_t>(dirty)]);
        }
    }

    // ------------------------------------------------------------ construction
    Solution construct(int K, int start_index, std::mt19937_64& rng) {
        Counts k{I.pm, I.dm};
        std::uniform_int_distribution<int> mode_dist(0, 2);
        const int qmode = start_index == 0 ? 0 : mode_dist(rng);
        std::vector<int> chosen;
        std::vector<long long> keys;
        keys.reserve(pool.comp_by_pair.size());
        for (const auto& entry : pool.comp_by_pair) keys.push_back(entry.first);
        std::sort(keys.begin(), keys.end());
        for (long long key : keys) {
            const int idx = select_comp(pool.comp_by_pair.at(key), qmode, rng, nullptr);
            if (idx < 0) continue;
            if (opt.construct_composite_filter == PalnsCompositeFilter::DeltaPositive &&
                pool.comps[static_cast<std::size_t>(idx)].delta <= kEps) {
                continue;
            }
            chosen.push_back(idx);
        }
        std::stable_sort(chosen.begin(), chosen.end(), [&](int a, int b) {
            const Composite& x = pool.comps[static_cast<std::size_t>(a)];
            const Composite& y = pool.comps[static_cast<std::size_t>(b)];
            if (std::fabs(x.delta - y.delta) > kEps) return x.delta > y.delta;
            return x.d < y.d;
        });
        std::vector<Trip> trips;
        for (int idx : chosen) {
            const Trip t = make_composite_trip(pool, idx);
            if (trip_valid(t, k)) {
                consume_trip(t, k);
                trips.push_back(t);
            }
        }
        std::vector<Trip> singles = assign_singles(I, k, 0, rng);
        trips.insert(trips.end(), singles.begin(), singles.end());

        PalnsConstructOrder order = opt.construct_order;
        if (order == PalnsConstructOrder::Mixed) {
            if (start_index == 0) order = PalnsConstructOrder::Lpt;
            else {
                const int m = mode_dist(rng);
                order = m == 0 ? PalnsConstructOrder::Lpt : (m == 1 ? PalnsConstructOrder::Random : PalnsConstructOrder::Regret);
            }
        }
        Solution s;
        s.routes.assign(static_cast<std::size_t>(K), Route{});
        if (order == PalnsConstructOrder::Regret) {
            insert_dynamic(s, trips, nullptr, true, false, false);
        } else {
            if (order == PalnsConstructOrder::Lpt) {
                std::stable_sort(trips.begin(), trips.end(), [](const Trip& a, const Trip& b) { return a.d > b.d; });
            } else {
                std::shuffle(trips.begin(), trips.end(), rng);
            }
            for (const Trip& t : trips) insert_trip(s, t, false, false);
        }
        refresh_totals(I, s);
        return s;
    }

    // ------------------------------------------------------------ destroy
    double removal_gain(const Solution& s, const TrajLoc& t, int traj) const {
        const Route& r = s.routes[static_cast<std::size_t>(t.route)];
        double d = 0.0, c = 0.0;
        int prev = -1;
        for (const Op& o : r.ops) {
            if (o.traj == traj) continue;
            const int l = I.loc(o);
            d += I.srv(o);
            if (prev >= 0) { d += I.tt(prev, l); c += I.cc(prev, l); }
            prev = l;
        }
        if (phase == 1) {
            return (excess(r.d, I.T) - excess(d, I.T)) + eta_cur * (r.d - d);
        }
        return (r.c - c) + lambda * (excess(r.d, I.T) - excess(d, I.T));
    }

    // ops: 0 random, 1 worst, 2 shaw, 3 extended associated, 4 overload
    std::unordered_set<int> destroy(const Solution& s, int op, int k, std::mt19937_64& rng, int& ref_class) {
        const std::unordered_map<int, TrajLoc> trajs = collect_trajs(s);
        std::vector<int> ids;
        ids.reserve(trajs.size());
        for (const auto& e : trajs) ids.push_back(e.first);
        std::sort(ids.begin(), ids.end());
        std::unordered_set<int> rem;
        if (ids.empty()) return rem;
        k = std::min<int>(k, static_cast<int>(ids.size()));
        std::uniform_real_distribution<double> U01(0.0, 1.0);
        auto pick_biased = [&](std::vector<int>& sorted, double power) {
            const double y = U01(rng);
            std::size_t idx = static_cast<std::size_t>(std::pow(y, power) * static_cast<double>(sorted.size()));
            if (idx >= sorted.size()) idx = sorted.size() - 1U;
            const int v = sorted[idx];
            sorted.erase(sorted.begin() + static_cast<std::ptrdiff_t>(idx));
            return v;
        };
        if (op == 4 && (phase == 2 || s.V <= kEps)) op = 1;
        if (op == 0) {
            std::shuffle(ids.begin(), ids.end(), rng);
            for (int i = 0; i < k; ++i) rem.insert(ids[static_cast<std::size_t>(i)]);
            ref_class = trajs.at(ids.front()).p;
        } else if (op == 1) {
            std::vector<std::pair<double, int>> g;
            for (int id : ids) g.push_back({removal_gain(s, trajs.at(id), id), id});
            std::sort(g.begin(), g.end(), [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second < b.second;
            });
            std::vector<int> sorted;
            for (const auto& x : g) sorted.push_back(x.second);
            for (int i = 0; i < k && !sorted.empty(); ++i) rem.insert(pick_biased(sorted, 3.0));
            ref_class = trajs.at(*rem.begin()).p;
        } else if (op == 2 || op == 4) {
            int seed = -1;
            if (op == 2) {
                std::uniform_int_distribution<std::size_t> dist(0, ids.size() - 1U);
                seed = ids[dist(rng)];
            } else {
                int worst_route = 0;
                double worst = -1.0;
                for (std::size_t r = 0; r < s.routes.size(); ++r) {
                    const double ex = excess(s.routes[r].d, I.T);
                    if (ex > worst) { worst = ex; worst_route = static_cast<int>(r); }
                }
                double best_gain = -kInf;
                for (int id : ids) {
                    const TrajLoc& t = trajs.at(id);
                    if (t.route != worst_route) continue;
                    const double g = removal_gain(s, t, id);
                    if (g > best_gain) { best_gain = g; seed = id; }
                }
                if (seed < 0) seed = ids.front();
            }
            rem.insert(seed);
            const int sp = trajs.at(seed).p;
            ref_class = sp;
            std::vector<std::pair<double, int>> g;
            for (int id : ids) {
                if (id == seed) continue;
                g.push_back({pool.sim[static_cast<std::size_t>(sp)][static_cast<std::size_t>(trajs.at(id).p)], id});
            }
            std::sort(g.begin(), g.end(), [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second < b.second;
            });
            std::vector<int> sorted;
            for (const auto& x : g) sorted.push_back(x.second);
            for (int i = 1; i < k && !sorted.empty(); ++i) rem.insert(pick_biased(sorted, 6.0));
        } else {
            std::shuffle(ids.begin(), ids.end(), rng);
            ref_class = trajs.at(ids.front()).p;
            std::vector<std::pair<int, int>> blocks;
            for (int i = 0; i < k; ++i) {
                const TrajLoc& t = trajs.at(ids[static_cast<std::size_t>(i)]);
                const Route& r = s.routes[static_cast<std::size_t>(t.route)];
                route_blocks(r.ops, blocks);
                for (const auto& b : blocks) {
                    if (t.pi >= b.first && t.pi < b.second) {
                        for (int x = b.first; x < b.second; ++x) rem.insert(r.ops[static_cast<std::size_t>(x)].traj);
                        break;
                    }
                }
            }
        }
        return rem;
    }

    // ------------------------------------------------------------ selection
    std::vector<Trip> select_candidates(const Counts& cnt, int qop, std::mt19937_64& rng) const {
        std::vector<Trip> out;
        std::vector<int> U;
        for (std::size_t p = 0; p < cnt.u.size(); ++p) {
            for (int m = 0; m < cnt.u[p]; ++m) U.push_back(static_cast<int>(p));
        }
        const int P = npc();
        for (std::size_t a = 0; a < U.size(); ++a) {
            for (std::size_t b = a + 1U; b < U.size(); ++b) {
                auto it = pool.comp_by_pair.find(pair_key(U[a], U[b], P));
                if (it == pool.comp_by_pair.end()) continue;
                const int idx = select_comp(it->second, qop, rng, &cnt.z);
                if (idx >= 0) out.push_back(make_composite_trip(pool, idx));
            }
        }
        std::mt19937_64& mrng = const_cast<std::mt19937_64&>(rng);
        std::vector<Trip> singles = assign_singles(I, cnt, qop, mrng, phase == 2 && v2.sel_cost);
        out.insert(out.end(), singles.begin(), singles.end());
        return out;
    }

    // ------------------------------------------------------------ repair
    // ops: 0 random, 1 shaw, 2 saving, 3 greedy, 4 regret-2
    void repair(Solution& s, std::vector<Trip> cand, Counts cnt, int rop, int ref_class, std::mt19937_64& rng) {
        const bool phase2 = phase == 2;
        if (rop <= 2) {
            std::vector<std::size_t> order(cand.size());
            std::iota(order.begin(), order.end(), 0U);
            if (rop == 0) {
                std::shuffle(order.begin(), order.end(), rng);
            } else if (rop == 1) {
                auto key = [&](const Trip& t) {
                    if (t.composite) return pool.sim[static_cast<std::size_t>(t.p0)][static_cast<std::size_t>(t.p1)];
                    const int r = ref_class >= 0 ? ref_class : t.p0;
                    return pool.sim[static_cast<std::size_t>(r)][static_cast<std::size_t>(t.p0)];
                };
                std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
                    return key(cand[a]) > key(cand[b]);
                });
            } else {
                std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
                    const Trip& x = cand[a];
                    const Trip& y = cand[b];
                    if (x.composite != y.composite) return x.composite;
                    if (x.composite) return x.delta > y.delta;
                    return x.d < y.d;
                });
            }
            for (std::size_t idx : order) {
                if (!trip_valid(cand[idx], cnt)) continue;
                insert_trip(s, cand[idx], true, phase2);
                consume_trip(cand[idx], cnt);
            }
        } else {
            insert_dynamic(s, cand, &cnt, rop == 4, true, phase2);
        }
        bool left = false;
        for (int u : cnt.u) if (u > 0) { left = true; break; }
        if (left) {
            std::vector<Trip> rest = assign_singles(I, cnt, 0, rng, phase == 2 && v2.sel_cost);
            std::stable_sort(rest.begin(), rest.end(), [](const Trip& a, const Trip& b) { return a.d > b.d; });
            for (const Trip& t : rest) insert_trip(s, t, true, phase2);
        }
        refresh_totals(I, s);
    }

    // ------------------------------------------------------------ local search
    // Objective change of replacing route values (d,c) by (nd,nc).
    double route_delta(double d, double c, double nd, double nc) const {
        if (phase == 1) {
            return (excess(nd, I.T) - excess(d, I.T)) + eta_cur * (nd - d);
        }
        return (nc - c) + lambda * (excess(nd, I.T) - excess(d, I.T));
    }

    void set_route(Solution& s, std::size_t r, std::vector<Op> ops) {
        Route& route = s.routes[r];
        const double od = route.d, oc = route.c;
        route.ops.swap(ops);
        refresh_route(I, route);
        s.V += excess(route.d, I.T) - excess(od, I.T);
        s.D += route.d - od;
        s.C += route.c - oc;
    }

    bool ls_dswap(Solution& s) {
        struct DPos { int r, i; };
        std::vector<std::vector<DPos>> by_type(static_cast<std::size_t>(I.ntype));
        for (std::size_t r = 0; r < s.routes.size(); ++r) {
            const auto& ops = s.routes[r].ops;
            for (std::size_t i = 0; i < ops.size(); ++i) {
                if (ops[i].kind == kDeliv) {
                    by_type[static_cast<std::size_t>(I.df[static_cast<std::size_t>(ops[i].cls)])].push_back(
                        DPos{static_cast<int>(r), static_cast<int>(i)});
                }
            }
        }
        double best = -1e-9;
        DPos ba{-1, -1}, bb{-1, -1};
        auto local = [&](const Route& route, int i, int nw, double& dd, double& dc) {
            const int m = static_cast<int>(route.ops.size());
            const int w = I.loc(route.ops[static_cast<std::size_t>(i)]);
            const int prev = i > 0 ? I.loc(route.ops[static_cast<std::size_t>(i - 1)]) : -1;
            const int next = i + 1 < m ? I.loc(route.ops[static_cast<std::size_t>(i + 1)]) : -1;
            dd = I.lt(prev, nw) + I.lt(nw, next) - I.lt(prev, w) - I.lt(w, next);
            dc = I.lc(prev, nw) + I.lc(nw, next) - I.lc(prev, w) - I.lc(w, next);
        };
        for (const auto& list : by_type) {
            for (std::size_t a = 0; a < list.size(); ++a) {
                for (std::size_t b = a + 1U; b < list.size(); ++b) {
                    const DPos x = list[a], y = list[b];
                    Route& rx = s.routes[static_cast<std::size_t>(x.r)];
                    Route& ry = s.routes[static_cast<std::size_t>(y.r)];
                    const int cx = rx.ops[static_cast<std::size_t>(x.i)].cls;
                    const int cy = ry.ops[static_cast<std::size_t>(y.i)].cls;
                    if (cx == cy) continue;
                    double delta;
                    if (x.r == y.r) {
                        std::vector<Op> ops = rx.ops;
                        std::swap(ops[static_cast<std::size_t>(x.i)].cls, ops[static_cast<std::size_t>(y.i)].cls);
                        double nd, nc;
                        eval_ops(I, ops, nd, nc);
                        delta = route_delta(rx.d, rx.c, nd, nc);
                    } else {
                        double d1, c1, d2, c2;
                        local(rx, x.i, I.dw[static_cast<std::size_t>(cy)], d1, c1);
                        local(ry, y.i, I.dw[static_cast<std::size_t>(cx)], d2, c2);
                        delta = route_delta(rx.d, rx.c, rx.d + d1, rx.c + c1) +
                                route_delta(ry.d, ry.c, ry.d + d2, ry.c + c2);
                    }
                    if (delta < best) {
                        best = delta; ba = x; bb = y;
                        if (opt.ls_first_improvement) goto apply_move;
                    }
                }
            }
        }
    apply_move:
        if (ba.r < 0) return false;
        {
            std::vector<Op> ox = s.routes[static_cast<std::size_t>(ba.r)].ops;
            if (ba.r == bb.r) {
                std::swap(ox[static_cast<std::size_t>(ba.i)].cls, ox[static_cast<std::size_t>(bb.i)].cls);
                set_route(s, static_cast<std::size_t>(ba.r), std::move(ox));
            } else {
                std::vector<Op> oy = s.routes[static_cast<std::size_t>(bb.r)].ops;
                std::swap(ox[static_cast<std::size_t>(ba.i)].cls, oy[static_cast<std::size_t>(bb.i)].cls);
                set_route(s, static_cast<std::size_t>(ba.r), std::move(ox));
                set_route(s, static_cast<std::size_t>(bb.r), std::move(oy));
            }
        }
        return true;
    }

    struct BlockRef {
        int r = -1, s = 0, e = 0;
        double d = 0.0, c = 0.0;
        int h = -1, g = -1;
    };

    void collect_blocks(const Solution& s, std::vector<std::vector<BlockRef>>& out) const {
        out.assign(s.routes.size(), {});
        std::vector<std::pair<int, int>> bl;
        for (std::size_t r = 0; r < s.routes.size(); ++r) {
            const auto& ops = s.routes[r].ops;
            route_blocks(ops, bl);
            for (const auto& b : bl) {
                BlockRef x;
                x.r = static_cast<int>(r); x.s = b.first; x.e = b.second;
                std::vector<Op> sub(ops.begin() + b.first, ops.begin() + b.second);
                eval_ops(I, sub, x.d, x.c);
                x.h = I.loc(ops[static_cast<std::size_t>(b.first)]);
                x.g = I.loc(ops[static_cast<std::size_t>(b.second - 1)]);
                out[r].push_back(x);
            }
        }
    }

    // Moves of whole blocks: inter-route swap and relocate, intra-route
    // relocate and block-order reversal.
    bool ls_blocks(Solution& s, int kind) {
        std::vector<std::vector<BlockRef>> bl;
        collect_blocks(s, bl);
        const std::size_t R = s.routes.size();
        double best = -1e-9;
        std::size_t br1 = 0, br2 = 0;
        std::vector<Op> bo1, bo2;
        bool found = false;
        auto locat = [&](std::size_t r, int i) {
            return (i < 0 || i >= static_cast<int>(s.routes[r].ops.size())) ? -1 : I.loc(s.routes[r].ops[static_cast<std::size_t>(i)]);
        };
        auto build_swap = [&](const BlockRef& a, const BlockRef& b, std::vector<Op>& out) {
            const auto& oa = s.routes[static_cast<std::size_t>(a.r)].ops;
            const auto& ob = s.routes[static_cast<std::size_t>(b.r)].ops;
            out.assign(oa.begin(), oa.begin() + a.s);
            out.insert(out.end(), ob.begin() + b.s, ob.begin() + b.e);
            out.insert(out.end(), oa.begin() + a.e, oa.end());
        };
        if (kind == 0) {
            for (std::size_t r1 = 0; r1 < R; ++r1) {
                for (std::size_t r2 = r1 + 1U; r2 < R; ++r2) {
                    const Route& A = s.routes[r1];
                    const Route& Bv = s.routes[r2];
                    for (const BlockRef& a : bl[r1]) {
                        const int pa = locat(r1, a.s - 1), na = locat(r1, a.e);
                        for (const BlockRef& b : bl[r2]) {
                            const int pb = locat(r2, b.s - 1), nb = locat(r2, b.e);
                            const double nd1 = A.d - a.d - I.lt(pa, a.h) - I.lt(a.g, na) + b.d + I.lt(pa, b.h) + I.lt(b.g, na);
                            const double nc1 = A.c - a.c - I.lc(pa, a.h) - I.lc(a.g, na) + b.c + I.lc(pa, b.h) + I.lc(b.g, na);
                            const double nd2 = Bv.d - b.d - I.lt(pb, b.h) - I.lt(b.g, nb) + a.d + I.lt(pb, a.h) + I.lt(a.g, nb);
                            const double nc2 = Bv.c - b.c - I.lc(pb, b.h) - I.lc(b.g, nb) + a.c + I.lc(pb, a.h) + I.lc(a.g, nb);
                            const double delta = route_delta(A.d, A.c, nd1, nc1) + route_delta(Bv.d, Bv.c, nd2, nc2);
                            if (delta < best) {
                                best = delta; found = true; br1 = r1; br2 = r2;
                                build_swap(a, b, bo1);
                                build_swap(b, a, bo2);
                                if (opt.ls_first_improvement) goto done;
                            }
                        }
                    }
                }
            }
        } else if (kind == 1) {
            for (std::size_t r1 = 0; r1 < R; ++r1) {
                const Route& A = s.routes[r1];
                for (const BlockRef& a : bl[r1]) {
                    const int pa = locat(r1, a.s - 1), na = locat(r1, a.e);
                    const double nd1 = A.d - a.d - I.lt(pa, a.h) - I.lt(a.g, na) + I.lt(pa, na);
                    const double nc1 = A.c - a.c - I.lc(pa, a.h) - I.lc(a.g, na) + I.lc(pa, na);
                    const double del1 = route_delta(A.d, A.c, nd1, nc1);
                    for (std::size_t r2 = 0; r2 < R; ++r2) {
                        if (r2 == r1) continue;
                        const Route& Bv = s.routes[r2];
                        std::vector<int> gaps{0};
                        for (const BlockRef& b : bl[r2]) gaps.push_back(b.e);
                        for (int g : gaps) {
                            const int pb = locat(r2, g - 1), nb = locat(r2, g);
                            const double nd2 = Bv.d + a.d + I.lt(pb, a.h) + I.lt(a.g, nb) - I.lt(pb, nb);
                            const double nc2 = Bv.c + a.c + I.lc(pb, a.h) + I.lc(a.g, nb) - I.lc(pb, nb);
                            const double delta = del1 + route_delta(Bv.d, Bv.c, nd2, nc2);
                            if (delta < best) {
                                best = delta; found = true; br1 = r1; br2 = r2;
                                const auto& oa = A.ops;
                                bo1.assign(oa.begin(), oa.begin() + a.s);
                                bo1.insert(bo1.end(), oa.begin() + a.e, oa.end());
                                const auto& ob = Bv.ops;
                                bo2.assign(ob.begin(), ob.begin() + g);
                                bo2.insert(bo2.end(), oa.begin() + a.s, oa.begin() + a.e);
                                bo2.insert(bo2.end(), ob.begin() + g, ob.end());
                                if (opt.ls_first_improvement) goto done;
                            }
                        }
                    }
                }
            }
        } else {
            for (std::size_t r = 0; r < R; ++r) {
                const Route& A = s.routes[r];
                const std::vector<BlockRef>& b = bl[r];
                const int nb = static_cast<int>(b.size());
                if (nb < 2) continue;
                auto assemble = [&](const std::vector<int>& order, std::vector<Op>& out) {
                    out.clear();
                    for (int idx : order) out.insert(out.end(), A.ops.begin() + b[static_cast<std::size_t>(idx)].s, A.ops.begin() + b[static_cast<std::size_t>(idx)].e);
                };
                std::vector<int> base(static_cast<std::size_t>(nb));
                std::iota(base.begin(), base.end(), 0);
                std::vector<Op> trial;
                // relocate block x to position y
                for (int x = 0; x < nb; ++x) {
                    for (int y = 0; y < nb; ++y) {
                        if (y == x) continue;
                        std::vector<int> order = base;
                        order.erase(order.begin() + x);
                        order.insert(order.begin() + y, x);
                        assemble(order, trial);
                        double nd, nc;
                        eval_ops(I, trial, nd, nc);
                        const double delta = route_delta(A.d, A.c, nd, nc);
                        if (delta < best) {
                            best = delta; found = true; br1 = r; br2 = r; bo1 = trial;
                            if (opt.ls_first_improvement) goto done;
                        }
                    }
                }
                // reverse a block-order segment
                for (int x = 0; x < nb; ++x) {
                    for (int y = x + 1; y < nb; ++y) {
                        std::vector<int> order = base;
                        std::reverse(order.begin() + x, order.begin() + y + 1);
                        assemble(order, trial);
                        double nd, nc;
                        eval_ops(I, trial, nd, nc);
                        const double delta = route_delta(A.d, A.c, nd, nc);
                        if (delta < best) {
                            best = delta; found = true; br1 = r; br2 = r; bo1 = trial;
                            if (opt.ls_first_improvement) goto done;
                        }
                    }
                }
            }
        }
    done:
        if (!found) return false;
        set_route(s, br1, bo1);
        if (br2 != br1) set_route(s, br2, bo2);
        return true;
    }

    // Removes one skip trajectory and reinserts it at its best position
    // (contiguous or split) over all routes. The route is modified in place
    // and restored after each evaluation.
    bool ls_relocate(Solution& s) {
        const std::unordered_map<int, TrajLoc> trajs = collect_trajs(s);
        std::vector<int> ids;
        ids.reserve(trajs.size());
        for (const auto& e : trajs) ids.push_back(e.first);
        std::sort(ids.begin(), ids.end());
        const double f0 = objective(s);
        double best_f = f0 - 1e-9;
        int best_id = -1;
        for (int id : ids) {
            const TrajLoc& t = trajs.at(id);
            const std::size_t ri = static_cast<std::size_t>(t.route);
            std::vector<Op> saved = s.routes[ri].ops;
            std::vector<Op> kept;
            kept.reserve(saved.size());
            Op op_p, op_e, op_d;
            for (const Op& o : saved) {
                if (o.traj != id) { kept.push_back(o); continue; }
                if (o.kind == kPick) op_p = o; else if (o.kind == kEmpty) op_e = o; else op_d = o;
            }
            set_route(s, ri, std::move(kept));
            BlockInfo B = make_block(I, std::vector<Op>{op_p, op_e, op_d});
            const InsEval e = eval_insertion(I, s, B, true, phase == 2, lambda, nullptr);
            double f = kInf;
            if (e.route >= 0) {
                const double V = s.V + e.dV, D = s.D + e.dd, C = s.C + e.dc;
                f = phase == 1 ? V + eta_cur * D : C + lambda * V;
            }
            if (f < best_f && opt.ls_first_improvement) {
                apply_insertion(I, s, B, e);
                refresh_totals(I, s);
                return true;
            }
            set_route(s, ri, std::move(saved));
            if (f < best_f) {
                best_f = f;
                best_id = id;
            }
        }
        if (best_id < 0) return false;
        const TrajLoc& t = trajs.at(best_id);
        const std::size_t ri = static_cast<std::size_t>(t.route);
        std::vector<Op> kept;
        Op op_p, op_e, op_d;
        for (const Op& o : s.routes[ri].ops) {
            if (o.traj != best_id) { kept.push_back(o); continue; }
            if (o.kind == kPick) op_p = o; else if (o.kind == kEmpty) op_e = o; else op_d = o;
        }
        set_route(s, ri, std::move(kept));
        BlockInfo B = make_block(I, std::vector<Op>{op_p, op_e, op_d});
        const InsEval e = eval_insertion(I, s, B, true, phase == 2, lambda, nullptr);
        apply_insertion(I, s, B, e);
        refresh_totals(I, s);
        return true;
    }

    // ---- version 2: trajectory relocate with a cache -----------------------
    // The same move as ls_relocate (same order, same tie-breaking, so the same
    // search), but the best insertion of a trajectory into a route is looked
    // up when that route has not changed since it was last computed:
    //   foreign route: depends on the route's operations (kind, class) and on
    //                  the classes of the trajectory;
    //   own route:     depends on the route's operations with their trajectory
    //                  labels and on the trajectory.
    // The cache is cleared whenever the objective changes (phase, lambda, eta).
    struct RelocEntry {
        InsEval e;
        double kd = 0.0, kc = 0.0;  // own route: duration and cost of the route without the trajectory
        int m = -1;                 // size of the route the entry was computed on (guard)
    };
    std::unordered_map<std::uint64_t, RelocEntry> reloc_cache;
    int reloc_phase = 0;
    double reloc_lambda = 0.0, reloc_eta = 0.0;
    long long stat_reloc_hit = 0, stat_reloc_miss = 0;
    std::vector<RouteCtx> reloc_ctx;       // contexts of the routes of the current scan
    std::vector<char> reloc_ctx_ok;
    RouteCtx reloc_own_ctx;
    SplitScratch reloc_scratch;
    SplitScratch ins_scratch;              // scratch of insert_dynamic

    static std::uint64_t route_hash(const Route& r, bool with_traj) {
        std::uint64_t h = 0x243F6A8885A308D3ULL;
        for (const Op& o : r.ops) {
            std::uint64_t lab = (static_cast<std::uint64_t>(o.kind + 1) << 40) ^
                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(o.cls)) << 8);
            if (with_traj) {
                lab ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(o.traj)) * 0x9E3779B97F4A7C15ULL;
            }
            h = splitmix(h ^ lab);
        }
        return h;
    }

    bool ls_relocate_cached(Solution& s) {
        if (reloc_phase != phase || reloc_lambda != lambda || reloc_eta != eta_cur ||
            reloc_cache.size() > 200000U) {
            reloc_cache.clear();
            reloc_phase = phase;
            reloc_lambda = lambda;
            reloc_eta = eta_cur;
        }
        const bool phase2 = phase == 2;
        const std::size_t R = s.routes.size();
        std::vector<std::uint64_t> hc(R), ht(R);
        reloc_ctx.resize(R);
        reloc_ctx_ok.assign(R, 0);
        for (std::size_t r = 0; r < R; ++r) {
            hc[r] = route_hash(s.routes[r], false);
            ht[r] = route_hash(s.routes[r], true);
        }
        const std::unordered_map<int, TrajLoc> trajs = collect_trajs(s);
        std::vector<int> ids;
        ids.reserve(trajs.size());
        for (const auto& e : trajs) ids.push_back(e.first);
        std::sort(ids.begin(), ids.end());
        const double f0 = objective(s);
        double best_f = f0 - 1e-9;
        int best_id = -1;
        InsEval best_e;
        auto without = [&](std::size_t ri, int id) {
            std::vector<Op> kept;
            kept.reserve(s.routes[ri].ops.size());
            for (const Op& o : s.routes[ri].ops) if (o.traj != id) kept.push_back(o);
            return kept;
        };
        auto apply = [&](int id, const InsEval& e) {
            const TrajLoc& t = trajs.at(id);
            const std::size_t ri = static_cast<std::size_t>(t.route);
            const Route& rr = s.routes[ri];
            BlockInfo B = make_block(I, std::vector<Op>{rr.ops[static_cast<std::size_t>(t.pi)],
                                                        rr.ops[static_cast<std::size_t>(t.ei)],
                                                        rr.ops[static_cast<std::size_t>(t.di)]});
            set_route(s, ri, without(ri, id));
            apply_insertion(I, s, B, e);
            refresh_totals(I, s);
        };
        for (int id : ids) {
            const TrajLoc& t = trajs.at(id);
            const std::size_t ri = static_cast<std::size_t>(t.route);
            const Route& rr = s.routes[ri];
            const BlockInfo B = make_block(I, std::vector<Op>{rr.ops[static_cast<std::size_t>(t.pi)],
                                                              rr.ops[static_cast<std::size_t>(t.ei)],
                                                              rr.ops[static_cast<std::size_t>(t.di)]});
            // own route without the trajectory
            const std::uint64_t okey = splitmix(ht[ri] ^ (0xD1B54A32D192ED03ULL *
                static_cast<std::uint64_t>(static_cast<std::uint32_t>(id) + 1U)));
            RelocEntry& own = reloc_cache[okey];
            const int own_m = static_cast<int>(rr.ops.size()) - 3;
            if (own.m != own_m) {
                Route kept;
                kept.ops = without(ri, id);
                refresh_route(I, kept);
                make_ctx(I, kept, reloc_own_ctx);
                own.e = eval_insertion_route(I, kept, reloc_own_ctx, static_cast<int>(ri), B, true, phase2,
                                             lambda, reloc_scratch);
                own.kd = kept.d;
                own.kc = kept.c;
                own.m = own_m;
                ++stat_reloc_miss;
            } else {
                ++stat_reloc_hit;
            }
            const double own_kd = own.kd, own_kc = own.kc;
            InsEval e;  // best over all routes, in route order as in eval_insertion
            const std::uint64_t cls_mix = 0x9E3779B97F4A7C15ULL * (static_cast<std::uint64_t>(static_cast<std::uint32_t>(t.p)) + 1U) ^
                0xC2B2AE3D27D4EB4FULL * (static_cast<std::uint64_t>(static_cast<std::uint32_t>(t.q)) + 1U);
            const InsEval own_e = own.e;  // copy: the map may rehash below
            for (std::size_t r = 0; r < R; ++r) {
                InsEval er;
                if (r == ri) {
                    er = own_e;
                } else {
                    const std::uint64_t fkey = splitmix(hc[r] ^ cls_mix) ^ 0x5851F42D4C957F2DULL;
                    RelocEntry& en = reloc_cache[fkey];
                    const int fm = static_cast<int>(s.routes[r].ops.size());
                    if (en.m != fm) {
                        if (!reloc_ctx_ok[r]) {
                            make_ctx(I, s.routes[r], reloc_ctx[r]);
                            reloc_ctx_ok[r] = 1;
                        }
                        en.e = eval_insertion_route(I, s.routes[r], reloc_ctx[r], static_cast<int>(r), B, true,
                                                    phase2, lambda, reloc_scratch);
                        en.m = fm;
                        ++stat_reloc_miss;
                    } else {
                        ++stat_reloc_hit;
                    }
                    er = en.e;
                }
                er.route = er.route < 0 ? -1 : static_cast<int>(r);
                if (ins_better(er, e, phase2, lambda)) e = er;
            }
            double f = kInf;
            if (e.route >= 0) {
                // totals after the removal, computed as set_route does
                const double V0 = s.V + (excess(own_kd, I.T) - excess(rr.d, I.T));
                const double D0 = s.D + (own_kd - rr.d);
                const double C0 = s.C + (own_kc - rr.c);
                const double V = V0 + e.dV, D = D0 + e.dd, C = C0 + e.dc;
                f = phase == 1 ? V + eta_cur * D : C + lambda * V;
            }
            if (f < best_f && opt.ls_first_improvement) {
                apply(id, e);
                return true;
            }
            if (f < best_f) {
                best_f = f;
                best_id = id;
                best_e = e;
            }
        }
        if (best_id < 0) return false;
        apply(best_id, best_e);
        return true;
    }

    // Version 2 profile (config profile=1): seconds and calls per part of an iteration.
    struct Profile {
        double ruin = 0.0, dswap = 0.0, blocks = 0.0, reloc = 0.0;
        long long n_ruin = 0, n_dswap = 0, n_blocks = 0, n_reloc = 0, ls_calls = 0, ls_moves = 0;
    };
    Profile prof;

    // Version 2 diagnostics (config diag=1): Phase II counters. Nothing here is
    // read by the search.
    struct Diag {
        // per operator: uses, accepted candidates, candidates better than the
        // current solution, new thread bests
        long long d_use[5] = {}, d_acc[5] = {}, d_imp[5] = {}, d_best[5] = {};
        long long q_use[3] = {}, q_acc[3] = {}, q_imp[3] = {}, q_best[3] = {};
        long long r_use[5] = {}, r_acc[5] = {}, r_imp[5] = {}, r_best[5] = {};
        // removal size (index min(k, 15)): uses and new thread bests
        long long k_use[16] = {}, k_best[16] = {};
        long long iters = 0, feasible = 0, seen_before = 0, accepted = 0, uphill = 0;
        // repair: iterations with a composite candidate, composite candidates,
        // trips inserted (composite, single, single in split form)
        long long comp_iters = 0, comp_cand = 0, comp_ins = 0, single_ins = 0, single_split = 0;
        // local search moves by type, and the stage after which a new thread
        // best was first reached
        long long ls_dswap = 0, ls_blocks[3] = {}, ls_reloc = 0;
        long long best_after_repair = 0, best_after_light = 0, best_after_deep = 0;
        long long lam_back = 0;  // returns to a feasible solution (option lam_back)
        long long indep_adopt = 0;  // shared best taken by an independent thread (option indep_adopt)
        struct Cycle {
            double start = 0.0, low = kInf, gap_sum = 0.0, lam = 0.0;
            long long n = 0, acc = 0, feas = 0, best = 0;
            // current solution: iterations in which it is infeasible, their
            // total excess, and the range of the penalty weight
            long long cur_inf = 0;
            double cur_inf_V = 0.0, lam_lo = kInf, lam_hi = 0.0;
            double mul = 1.0;  // temperature multiplier of the cycle (option reheat_max)
        };
        std::vector<Cycle> cycles;
    };
    Diag diag;

    // Version 2 diagnostics (config pool_ratio, pool_file): blocks of the good
    // feasible Phase II candidates of this thread. Nothing here is read by the
    // search. A block is keyed by its operations (kind, class), 3 bytes each.
    struct PoolEvent {
        int block = 0;
        long long gi = 0;
        double cost = 0.0;
    };
    std::unordered_map<std::string, int> pool_index;
    std::vector<std::string> pool_keys;
    std::vector<double> pool_min;        // cheapest solution (original cost) containing the block
    std::vector<PoolEvent> pool_events;  // block first seen, or its cheapest solution got cheaper (only with pool_file)
    long long pool_solutions = 0;
    // Pool of the last solve_fixed_k call, merged over the threads (master object only).
    std::vector<PoolBlock> last_pool;

    void pool_note(const Solution& s) {
        const double c = original_cost(s);
        ++pool_solutions;
        std::string key;
        for (const Route& r : s.routes) {
            int load = 0;
            key.clear();
            for (const Op& o : r.ops) {
                key.push_back(static_cast<char>(o.kind));
                key.push_back(static_cast<char>(o.cls & 0xFF));
                key.push_back(static_cast<char>((o.cls >> 8) & 0xFF));
                load += load_delta(o);
                if (load != 0) continue;
                auto it = pool_index.find(key);
                if (it == pool_index.end()) {
                    it = pool_index.emplace(key, static_cast<int>(pool_keys.size())).first;
                    pool_keys.push_back(key);
                    pool_min.push_back(kInf);
                }
                const std::size_t id = static_cast<std::size_t>(it->second);
                if (c < pool_min[id] - 1e-9) {
                    pool_min[id] = c;
                    if (!v2.pool_file.empty()) pool_events.push_back({it->second, gi + 1, c});
                }
                key.clear();
            }
        }
    }

    void local_search(Solution& s, bool deep) {
        int moves = 0;
        bool improved = true;
        if (v2.profile) {
            ++prof.ls_calls;
            auto timed = [&](double& acc, long long& n, bool r, Clock::time_point t) {
                acc += std::chrono::duration<double>(Clock::now() - t).count();
                ++n;
                return r;
            };
            while (improved && moves < 200 && !timed_out()) {
                improved = false;
                auto t = Clock::now();
                bool hit = timed(prof.dswap, prof.n_dswap, ls_dswap(s), t);
                for (int kind = 0; kind < 3 && !hit; ++kind) {
                    t = Clock::now();
                    hit = timed(prof.blocks, prof.n_blocks, ls_blocks(s, kind), t);
                }
                if (!hit && deep) {
                    t = Clock::now();
                    hit = timed(prof.reloc, prof.n_reloc,
                                v2.fast_relocate ? ls_relocate_cached(s) : ls_relocate(s), t);
                }
                if (hit) {
                    improved = true;
                    ++moves;
                    ++prof.ls_moves;
                }
            }
            refresh_totals(I, s);
            return;
        }
        if (v2.diag && phase == 2) {
            // same moves in the same order as below, counted by type
            while (improved && moves < 200 && !timed_out()) {
                improved = false;
                bool hit = ls_dswap(s);
                if (hit) ++diag.ls_dswap;
                for (int kind = 0; kind < 3 && !hit; ++kind) {
                    hit = ls_blocks(s, kind);
                    if (hit) ++diag.ls_blocks[kind];
                }
                if (!hit && deep) {
                    hit = v2.fast_relocate ? ls_relocate_cached(s) : ls_relocate(s);
                    if (hit) ++diag.ls_reloc;
                }
                if (hit) {
                    improved = true;
                    ++moves;
                }
            }
            refresh_totals(I, s);
            return;
        }
        while (improved && moves < 200 && !timed_out()) {
            improved = false;
            if (ls_dswap(s) || ls_blocks(s, 0) || ls_blocks(s, 1) || ls_blocks(s, 2) ||
                (deep && (v2.fast_relocate ? ls_relocate_cached(s) : ls_relocate(s)))) {
                improved = true;
                ++moves;
            }
        }
        refresh_totals(I, s);
    }

    // ------------------------------------------------------------ search
    struct SearchState {
        bool have_feasible = false;
        Solution best_feasible;
        double best_violation = kInf;
        double time_to_feasible = -1.0;
        std::vector<std::pair<double, double>> trace;  // (time, original cost)
        // Least-violated solution, lexicographic (V, C), kept only while no
        // feasible solution exists (diagnostics; does not affect the search).
        bool have_least = false;
        Solution least;
    };

    void record(SearchState& st, const Solution& s) {
        st.best_violation = std::min(st.best_violation, s.V);
        if (!st.have_feasible &&
            (!st.have_least || s.V < st.least.V - kEps ||
             (std::fabs(s.V - st.least.V) <= kEps && s.C < st.least.C - kEps))) {
            st.least = s;
            st.have_least = true;
        }
        if (s.V > kEps) return;
        publish(s);
        if (!st.have_feasible || s.C < st.best_feasible.C - kEps) {
            st.best_feasible = s;
            st.have_feasible = true;
            st.trace.push_back({elapsed(), original_cost(s)});
        }
        if (st.time_to_feasible < 0.0) st.time_to_feasible = elapsed();
    }

    long long run_phase(int which, Solution start, SearchState& st, std::mt19937_64& rng) {
        // (definition below)
        return run_phase_impl(which, std::move(start), st, rng);
    }

    // Rescue neighbourhood. With probability 0.2 a single random route (its
    // full rebuild lowers the total duration). Otherwise an overloaded seed
    // route and up to rescue_routes - 1 routes with slack (d < T), chosen
    // with a bias towards the closest service locations.
    std::vector<std::size_t> rescue_neighbourhood(const Solution& s, std::mt19937_64& rng) const {
        std::uniform_real_distribution<double> U01(0.0, 1.0);
        const std::size_t R = s.routes.size();
        std::vector<std::size_t> over;
        for (std::size_t r = 0; r < R; ++r) if (s.routes[r].d > I.T + kEps) over.push_back(r);
        std::size_t seed;
        if (over.empty() || U01(rng) < 0.2) {
            std::uniform_int_distribution<std::size_t> d(0, R - 1U);
            return {d(rng)};
        }
        {
            std::uniform_int_distribution<std::size_t> d(0, over.size() - 1U);
            seed = over[d(rng)];
        }
        const int m = std::max(1, std::min<int>(opt.rescue_routes, static_cast<int>(R)));
        std::vector<int> ls;
        for (const Op& o : s.routes[seed].ops) ls.push_back(I.loc(o));
        std::vector<std::pair<double, std::size_t>> cand;
        for (std::size_t r = 0; r < R; ++r) {
            if (r == seed || s.routes[r].d >= I.T - kEps) continue;
            double best = kInf;
            for (const Op& o : s.routes[r].ops) {
                const int l = I.loc(o);
                for (int x : ls) best = std::min(best, std::min(I.tt(x, l), I.tt(l, x)));
            }
            cand.push_back({best, r});
        }
        std::sort(cand.begin(), cand.end());
        std::vector<std::size_t> out{seed};
        while (static_cast<int>(out.size()) < m && !cand.empty()) {
            const double y = U01(rng);
            std::size_t idx = static_cast<std::size_t>(y * y * static_cast<double>(cand.size()));
            if (idx >= cand.size()) idx = cand.size() - 1U;
            out.push_back(cand[idx].second);
            cand.erase(cand.begin() + static_cast<std::ptrdiff_t>(idx));
        }
        return out;
    }

    // Focus neighbourhood (version 2): a seed route (the routes take turns; ties
    // are broken at random) and focus_routes - 1 other routes, chosen with a
    // bias towards routes whose service locations are close to those of the
    // seed route (mean nearest travel time).
    std::vector<std::size_t> focus_neighbourhood(const Solution& s, std::mt19937_64& rng) {
        std::uniform_real_distribution<double> U01(0.0, 1.0);
        const std::size_t R = s.routes.size();
        if (focus_last.size() != R) focus_last.assign(R, -1);
        long long oldest = std::numeric_limits<long long>::max();
        for (long long v : focus_last) oldest = std::min(oldest, v);
        std::vector<std::size_t> ties;
        for (std::size_t r = 0; r < R; ++r) if (focus_last[r] == oldest) ties.push_back(r);
        std::uniform_int_distribution<std::size_t> d(0, ties.size() - 1U);
        const std::size_t seed = ties[d(rng)];
        focus_last[seed] = focus_count++;
        const int m = std::max(1, std::min<int>(v2.focus_routes, static_cast<int>(R)));
        std::vector<int> ls;
        for (const Op& o : s.routes[seed].ops) ls.push_back(I.loc(o));
        std::vector<std::pair<double, std::size_t>> cand;
        for (std::size_t r = 0; r < R; ++r) {
            if (r == seed) continue;
            double sum = 0.0;
            for (int x : ls) {
                double best = kInf;
                for (const Op& o : s.routes[r].ops) {
                    const int l = I.loc(o);
                    best = std::min(best, std::min(I.tt(x, l), I.tt(l, x)));
                }
                if (best < kInf) sum += best;
            }
            cand.push_back({ls.empty() ? 0.0 : sum / static_cast<double>(ls.size()), r});
        }
        std::sort(cand.begin(), cand.end());
        std::vector<std::size_t> out{seed};
        while (static_cast<int>(out.size()) < m && !cand.empty()) {
            const double y = U01(rng);
            std::size_t idx = static_cast<std::size_t>(y * y * static_cast<double>(cand.size()));
            if (idx >= cand.size()) idx = cand.size() - 1U;
            out.push_back(cand[idx].second);
            cand.erase(cand.begin() + static_cast<std::ptrdiff_t>(idx));
        }
        return out;
    }

    // Starts a focus episode on the thread solution S (feasible): at most
    // max_len iterations on focus_routes routes, all other routes fixed.
    void start_episode(const Solution& S, long long max_len, double start_ratio, std::mt19937_64& rng) {
        ep.nb = focus_neighbourhood(S, rng);
        ep.cur.routes.clear();
        for (std::size_t r : ep.nb) ep.cur.routes.push_back(S.routes[r]);
        refresh_totals(I, ep.cur);
        ep.best = ep.cur;
        ep.fcur = objective(ep.cur);
        ep.theta0 = std::max(1e-9, start_ratio * std::max(ep.fcur, 1e-9) / std::log(2.0));
        ep.len = std::max<long long>(1, std::min<long long>(std::max(1, v2.focus_len), max_len));
        ep.done = 0;
        ep.active = true;
        ++stat_focus_episodes;
    }

    // Ends the episode: its best sub-solution replaces the routes of S.
    void finish_episode(Solution& S, Solution& best, double& fS, double& fBest, SearchState& st) {
        for (std::size_t i = 0; i < ep.nb.size(); ++i) S.routes[ep.nb[i]] = ep.best.routes[i];
        refresh_totals(I, S);
        fS = objective(S);
        if (S.V <= kEps && S.C < fBest - 1e-9) {
            best = S;
            fBest = S.C;
            record(st, S);
            ++stat_focus_best;
        }
        ep.active = false;
    }

    long long run_phase_impl(int which, Solution start, SearchState& st, std::mt19937_64& rng) {
        phase = which;
        eta_cur = opt.eta;
        const bool oscillate = which == 1 && opt.phase1_eta_max > opt.eta;
        refresh_totals(I, start);
        if (which == 2) {
            lambda = std::max(1.0, 2.0 * start.C / std::max(1.0, start.D));
        }
        Solution S = start;
        Solution best = start;
        double fS = objective(S);
        double fBest = which == 1 ? fS : start.C;
        const bool pool_on = which == 2 && (v2.pool_ratio > 0.0 || v2.pool_delta > 0.0);
        if (pool_on && start.V <= kEps) pool_note(start);
        int lam_window = 0, lam_feasible = 0;
        const double lam_min = v2.lam_lo > 0.0 ? v2.lam_lo * lambda : 1e-3;
        const double lam_max = std::max(lam_min, v2.lam_hi > 0.0 ? v2.lam_hi * lambda : 1e6);
        const double lam_step = v2.lam_step > 1.0 ? v2.lam_step : 1.5;
        const double lam_dstep = v2.lam_dstep > 1.0 ? v2.lam_dstep : lam_step;
        const int lam_win = std::max(1, v2.lam_win);
        long long inf_run = 0;  // consecutive iterations with an infeasible current solution
        const bool reheat = which == 2 && v2.reheat_max > 1.0;
        double theta_mul = 1.0;       // temperature multiplier of the current cycle
        double cycle_best = fBest;    // thread best at the start of the current cycle
        int stale_cycles = 0;
        int own_stale = 0;            // cycles in a row without a new best of this thread
        const bool keep_feas = which == 2 && v2.lam_back > 0 && v2.lam_to == 1;
        Solution Sfeas;         // last feasible current solution (lam_to=1)
        bool has_feas = false;
        if (which == 2) lambda = std::min(lam_max, std::max(lam_min, lambda));
        std::vector<double> wd(5, 1.0), wq(3, 1.0), wr(5, 1.0);
        std::vector<double> sd(5, 0.0), sq(3, 0.0), sr(5, 0.0);
        std::vector<int> nd(5, 0), nq(3, 0), nr(5, 0);
        const int weight_period = 50;
        const double rho = 0.1;
        const double start_ratio = (which == 1 && opt.sa_start_ratio_phase1 >= 0.0)
            ? opt.sa_start_ratio_phase1 : opt.sa_start_ratio;
        double theta0 = std::max(1e-9, start_ratio * std::max(fS, 1e-9) / std::log(2.0));
        const long long cycle = std::max<long long>(1,
            (which == 2 && opt.sa_cycle_iters_phase2 > 0) ? opt.sa_cycle_iters_phase2 : opt.sa_cycle_iters);
        std::unordered_set<std::uint64_t> seen;
        seen.insert(solution_hash(S));
        const int n = I.n;
        const int kmax = std::max(1, std::min(n, std::max(opt.kmax_min,
            static_cast<int>(std::ceil(opt.kmax_ratio * static_cast<double>(n))))));
        std::uniform_int_distribution<int> kdist(1, kmax);
        std::uniform_real_distribution<double> U01(0.0, 1.0);
        long long iter = 0;
        bool rescue = false;
        const bool focus_enabled = which == 2 && v2.focus_period > 0 && S.routes.size() >= 2U;
        ep.active = false;
        while (!timed_out()) {
            if (which == 1 && st.have_feasible) break;
            ++iter;
            if (oscillate && !rescue) {
                const double pos = static_cast<double>((iter - 1) % cycle) / static_cast<double>(cycle);
                eta_cur = std::exp(std::log(opt.phase1_eta_max) + (std::log(opt.eta) - std::log(opt.phase1_eta_max)) * pos);
                fS = objective(S);
                fBest = objective(best);
                if ((iter - 1) % cycle == 0) {
                    theta0 = std::max(1e-9, start_ratio * std::max(fS, 1e-9) / std::log(2.0));
                }
            }
            std::vector<double> wdd = wd;
            if (which == 2) wdd[4] = 0.0;
            const int dop = roulette_index(wdd, rng);
            const int qop = roulette_index(wq, rng);
            const int rop = roulette_index(wr, rng);
            int k = kdist(rng);
            double reward = 0.0;
            double cand_V = 0.0;  // duration excess of the candidate of this iteration
            // Version 2: every focus_period-th SA cycle of Phase II is a focus cycle.
            const bool focus_cycle = focus_enabled && !rescue &&
                (((iter - 1) / cycle) % v2.focus_period) == v2.focus_period - 1;
            if (focus_cycle) {
                if (!ep.active) {
                    if (S.V > kEps) { S = best; fS = objective(S); }
                    start_episode(S, cycle - ((iter - 1) % cycle),
                                  v2.focus_start_ratio >= 0.0 ? v2.focus_start_ratio : start_ratio, rng);
                }
                Solution sps = ep.cur;
                {
                    std::size_t nops = 0;
                    for (const Route& r : sps.routes) nops += r.ops.size();
                    const int nsub = static_cast<int>(nops / 3U);
                    const int ksub = std::max(1, std::max(opt.kmax_min,
                        static_cast<int>(std::ceil(v2.focus_kmax_ratio * static_cast<double>(nsub)))));
                    std::uniform_int_distribution<int> kd(1, ksub);
                    k = kd(rng);
                    if (U01(rng) < v2.focus_rebuild) k = std::numeric_limits<int>::max();
                }
                int ref_class = -1;
                const std::unordered_set<int> rem = destroy(sps, dop, k, rng, ref_class);
                Counts cnt{std::vector<int>(static_cast<std::size_t>(npc()), 0),
                           std::vector<int>(static_cast<std::size_t>(ndc()), 0)};
                remove_trajs(I, sps, rem, cnt);
                std::vector<Trip> cand = select_candidates(cnt, qop, rng);
                repair(sps, std::move(cand), cnt, rop, ref_class, rng);
                local_search(sps, false);
                double fs = objective(sps);
                const double fb = objective(ep.best);
                if (v2.focus_deep == 1 ||
                    (v2.focus_deep == 2 && fs <= fb + opt.ls_deep_ratio * std::fabs(fb) + 1e-9)) {
                    local_search(sps, true);
                    fs = objective(sps);
                }
                cand_V = sps.V;
                const double delta = fs - ep.fcur;
                bool accept = delta < -1e-12;
                if (!accept) {
                    const double pos = static_cast<double>(ep.done) / static_cast<double>(ep.len);
                    const double theta = ep.theta0 * std::pow(opt.sa_end_ratio, pos);
                    accept = U01(rng) < std::exp(-delta / std::max(theta, 1e-300));
                }
                if (sps.V <= kEps && sps.C < ep.best.C - 1e-9) {
                    ep.best = sps;
                    reward = 10.0;
                } else if (delta < -1e-9) {
                    reward = 3.0;
                }
                if (accept) {
                    ep.cur = std::move(sps);
                    ep.fcur = fs;
                }
                ++ep.done;
                ++stat_focus_iters;
                if (ep.done >= ep.len) finish_episode(S, best, fS, fBest, st);
            } else {
                Solution sp = S;
                int ref_class = -1;
                // Rescue: operate on a sub-solution of a few routes, the others fixed.
                std::vector<std::size_t> nb;
                Solution sub;
                if (rescue) {
                    // keep a neighbourhood for 100 iterations (intensification)
                    if (rescue_nb_left <= 0 || rescue_nb.empty()) {
                        rescue_nb = rescue_neighbourhood(S, rng);
                        rescue_nb_left = 100;
                    }
                    --rescue_nb_left;
                    nb = rescue_nb;
                    for (std::size_t r : nb) sub.routes.push_back(S.routes[r]);
                    refresh_totals(I, sub);
                    if (U01(rng) < 1.0 / 3.0) k = std::numeric_limits<int>::max();  // rebuild all
                }
                Solution& tgt = rescue ? sub : sp;
                const auto ruin_t0 = v2.profile ? Clock::now() : Clock::time_point();
                const std::unordered_set<int> rem = destroy(tgt, dop, k, rng, ref_class);
                Counts cnt{std::vector<int>(static_cast<std::size_t>(npc()), 0),
                           std::vector<int>(static_cast<std::size_t>(ndc()), 0)};
                remove_trajs(I, tgt, rem, cnt);
                std::vector<Trip> cand = select_candidates(cnt, qop, rng);
                const bool dg = v2.diag && which == 2 && !rescue;
                if (dg) {
                    long long nc = 0;
                    for (const Trip& t : cand) if (t.composite) ++nc;
                    diag.comp_cand += nc;
                    if (nc > 0) ++diag.comp_iters;
                }
                repair(tgt, std::move(cand), cnt, rop, ref_class, rng);
                const bool dg_rep = dg && sp.V <= kEps && sp.C < fBest - 1e-9;
                bool dg_light = false;
                if (v2.profile) {
                    prof.ruin += std::chrono::duration<double>(Clock::now() - ruin_t0).count();
                    ++prof.n_ruin;
                }
                if (rescue) {
                    local_search(sub, false);
                    for (std::size_t i = 0; i < nb.size(); ++i) sp.routes[nb[i]] = sub.routes[i];
                    refresh_totals(I, sp);
                }
                bool run_ls = true;
                if (opt.ls_frequency == PalnsLsFrequency::OnImprove) run_ls = objective(sp) < fS - 1e-9;
                else if (opt.ls_frequency == PalnsLsFrequency::EveryM) run_ls = (iter % std::max(1, opt.ls_every_m)) == 0;
                if (run_ls) {
                    local_search(sp, false);
                    const double fl = objective(sp);
                    if (dg) dg_light = sp.V <= kEps && sp.C < fBest - 1e-9;
                    const double ref = which == 1 ? fBest : fBest;
                    if (fl <= ref + opt.ls_deep_ratio * std::fabs(ref) + 1e-9) local_search(sp, true);
                }
                refresh_totals(I, sp);
                record(st, sp);
                cand_V = sp.V;
                if (pool_on && sp.V <= kEps &&
                    sp.C <= std::max(fBest * (1.0 + v2.pool_ratio), fBest + v2.pool_delta) + 1e-9) pool_note(sp);
                const double fsp = objective(sp);
                {
                    const double delta = fsp - fS;
                    bool accept = false;
                    if (delta < -1e-12) {
                        accept = true;
                    } else {
                        const double pos = static_cast<double>((iter - 1) % cycle) / static_cast<double>(cycle);
                        const double theta = theta0 * theta_mul * std::pow(opt.sa_end_ratio, pos);
                        accept = U01(rng) < std::exp(-delta / std::max(theta, 1e-300));
                    }
                    const std::uint64_t h = solution_hash(sp);
                    const bool is_new = seen.find(h) == seen.end();
                    const bool new_best = which == 1
                        ? (oscillate ? (sp.V < best.V - kEps || (std::fabs(sp.V - best.V) <= kEps && sp.D < best.D - kEps))
                                     : (fsp < fBest - 1e-9))
                        : (sp.V <= kEps && sp.C < fBest - 1e-9);
                    if (new_best) {
                        best = sp;
                        fBest = which == 1 ? fsp : sp.C;
                        reward = 10.0;
                        if (which == 2) {
                            publish(sp);
                            ++stat_global_best;
                        }
                        if (rescue) { rbest = best; rbest_has = true; }
                    } else {
                        if (fsp < fS - 1e-9) reward = 3.0;
                        else if (accept && is_new) reward = 1.0;
                    }
                    if (dg) {
                        const std::size_t di = static_cast<std::size_t>(dop), qi = static_cast<std::size_t>(qop),
                                          ri = static_cast<std::size_t>(rop);
                        const std::size_t ki = static_cast<std::size_t>(std::min(k, 15));
                        const bool imp = fsp < fS - 1e-9;
                        ++diag.iters; ++diag.d_use[di]; ++diag.q_use[qi]; ++diag.r_use[ri]; ++diag.k_use[ki];
                        if (sp.V <= kEps) ++diag.feasible;
                        if (!is_new) ++diag.seen_before;
                        if (accept) { ++diag.accepted; ++diag.d_acc[di]; ++diag.q_acc[qi]; ++diag.r_acc[ri]; }
                        if (accept && delta > 1e-12) ++diag.uphill;
                        if (imp) { ++diag.d_imp[di]; ++diag.q_imp[qi]; ++diag.r_imp[ri]; }
                        if (new_best) {
                            ++diag.d_best[di]; ++diag.q_best[qi]; ++diag.r_best[ri]; ++diag.k_best[ki];
                            if (dg_rep) ++diag.best_after_repair;
                            else if (dg_light) ++diag.best_after_light;
                            else ++diag.best_after_deep;
                        }
                        const std::size_t ci = static_cast<std::size_t>((iter - 1) / cycle);
                        if (diag.cycles.size() <= ci) {
                            diag.cycles.resize(ci + 1U);
                            diag.cycles[ci].start = fBest;
                        }
                        Diag::Cycle& cy = diag.cycles[ci];
                        ++cy.n;
                        if (sp.V <= kEps) { ++cy.feas; cy.low = std::min(cy.low, sp.C); }
                        if (accept) ++cy.acc;
                        if (new_best) ++cy.best;
                        cy.gap_sum += ((accept ? fsp : fS) - fBest) / std::max(1.0, fBest);
                        cy.lam = lambda;
                        cy.mul = theta_mul;
                        cy.lam_lo = std::min(cy.lam_lo, lambda);
                        cy.lam_hi = std::max(cy.lam_hi, lambda);
                        const double cur_V = accept ? sp.V : S.V;
                        if (cur_V > kEps) { ++cy.cur_inf; cy.cur_inf_V += cur_V; }
                    }
                    if (accept) {
                        if (keep_feas && S.V <= kEps && sp.V > kEps) {
                            Sfeas = S;
                            has_feas = true;
                        }
                        S = sp;
                        fS = fsp;
                        seen.insert(h);
                    }
                }
            }
            if (which == 2) {
                ++lam_window;
                const double cur_V = ep.active ? ep.cur.V : S.V;
                if ((v2.lam_sig == 1 ? cur_V : cand_V) <= kEps) ++lam_feasible;
                if (lam_window >= lam_win) {
                    const double frac = static_cast<double>(lam_feasible) / lam_window;
                    if (frac < v2.lam_tlo) lambda = std::min(lam_max, lambda * lam_step);
                    else if (frac > v2.lam_thi) lambda = std::max(lam_min, lambda / lam_dstep);
                    lam_window = 0;
                    lam_feasible = 0;
                    fS = objective(S);
                    if (ep.active) ep.fcur = objective(ep.cur);
                }
                if (v2.lam_back > 0 && !ep.active) {
                    inf_run = S.V > kEps ? inf_run + 1 : 0;
                    if (inf_run >= v2.lam_back) {
                        S = (keep_feas && has_feas) ? Sfeas : best;
                        fS = objective(S);
                        inf_run = 0;
                        if (v2.diag) ++diag.lam_back;
                    }
                }
            }
            const bool cycle_end = iter % cycle == 0;
            if (which == 2 && cycle_end) own_stale = fBest < cycle_best - 1e-9 ? 0 : own_stale + 1;
            if (opt.sa_restart_best && cycle_end) {
                if (which == 2) {
                    Solution g;
                    const bool force = v2.indep_adopt > 0 && own_stale >= v2.indep_adopt;
                    if (adopt(g, fBest, force)) {
                        best = g;
                        fBest = g.C;
                        record(st, g);
                        own_stale = 0;
                        if (v2.diag && !adopt_shared) ++diag.indep_adopt;
                    }
                }
                S = best;
                fS = objective(S);
            }
            if (which == 2 && cycle_end) {
                if (reheat) {
                    if (fBest < cycle_best - 1e-9) {
                        theta_mul = 1.0;
                        stale_cycles = 0;
                    } else if (++stale_cycles >= std::max(1, v2.reheat_after)) {
                        theta_mul = std::min(v2.reheat_max, theta_mul + v2.reheat_step);
                        stale_cycles = 0;
                    }
                }
                cycle_best = fBest;
            }
            sd[static_cast<std::size_t>(dop)] += reward; ++nd[static_cast<std::size_t>(dop)];
            sq[static_cast<std::size_t>(qop)] += reward; ++nq[static_cast<std::size_t>(qop)];
            sr[static_cast<std::size_t>(rop)] += reward; ++nr[static_cast<std::size_t>(rop)];
            if (iter % weight_period == 0) {
                auto upd = [&](std::vector<double>& w, std::vector<double>& sc, std::vector<int>& cnts) {
                    for (std::size_t i = 0; i < w.size(); ++i) {
                        if (cnts[i] > 0) w[i] = (1.0 - rho) * w[i] + rho * (sc[i] / cnts[i]);
                        w[i] = std::max(w[i], 1e-3);
                        sc[i] = 0.0;
                        cnts[i] = 0;
                    }
                };
                upd(wd, sd, nd);
                upd(wq, sq, nq);
                upd(wr, sr, nr);
            }
            ++gi;
            if (hub != nullptr && gi % std::max(1, opt.sync_iters) == 0) {
                if (hub_arrive(*hub, wait_seconds)) break;
                // Phase I: even threads move on to Phase II with the shared best.
                if (which == 1) {
                    Solution g;
                    const bool force = hub->has &&
                        ((v2.p1_join >= 0 && gi - hub->first_i >= v2.p1_join) ||
                         (v2.rescue_p2 && hub->rescue));
                    if (adopt(g, kInf, force)) {
                        record(st, g);
                        break;
                    }
                    if (hub->rescue && !rescue) {
                        // Rescue: restart from this thread's least-violated solution.
                        rescue = true;
                        if (st.have_least) S = st.least;
                        eta_cur = opt.eta;  // rescue: the duration excess dominates
                        refresh_totals(I, S);
                        rescue_start_V = S.V;
                        best = S;
                        rbest = S;
                        rbest_has = true;
                        fS = objective(S);
                        fBest = fS;
                        theta0 = std::max(1e-9, start_ratio * std::max(fS, 1e-9) / std::log(2.0));
                        ++rescue_switches;
                    } else if (rescue && adopt_shared && hub->rhas &&
                               (hub->rsol.V < best.V - kEps ||
                                (std::fabs(hub->rsol.V - best.V) <= kEps && hub->rsol.D < best.D - kEps))) {
                        // Rescue: even threads continue from the best infeasible solution.
                        S = hub->rsol;
                        renumber(S);
                        refresh_totals(I, S);
                        best = S;
                        rbest = S;
                        fS = objective(S);
                        fBest = fS;
                        rescue_nb_left = 0;
                    }
                }
            }
        }
        if (ep.active) finish_episode(S, best, fS, fBest, st);  // stopped in the middle of an episode
        return iter;
    }

    // ------------------------------------------------------------ output
    bool build_action_routes(Solution s, int K, std::vector<CpActionRoute>& out, std::string& err) const {
        std::vector<std::vector<Op>> rs;
        for (const Route& r : s.routes) if (!r.ops.empty()) rs.push_back(r.ops);
        std::vector<std::pair<int, int>> bl;
        while (static_cast<int>(rs.size()) < K) {
            bool split = false;
            for (std::size_t r = 0; r < rs.size() && !split; ++r) {
                route_blocks(rs[r], bl);
                if (bl.size() >= 2U) {
                    std::vector<Op> tail(rs[r].begin() + bl.front().second, rs[r].end());
                    rs[r].resize(static_cast<std::size_t>(bl.front().second));
                    rs.push_back(std::move(tail));
                    split = true;
                }
            }
            if (!split) {
                err = "cannot fill the required number of nonempty routes";
                return false;
            }
        }
        if (static_cast<int>(rs.size()) > K) {
            err = "more nonempty routes than vehicles";
            return false;
        }
        // Zero-duration canonical swaps required by graph pruning rules
        // (41a: P->D, 41b: E->P, 41c: D->E at the same location).
        for (std::vector<Op>& ops : rs) {
            for (int guard = 0; guard < 100000; ++guard) {
                bool changed = false;
                for (std::size_t i = 0; i + 1U < ops.size(); ++i) {
                    const Op& a = ops[i];
                    const Op& b = ops[i + 1U];
                    if (I.loc(a) != I.loc(b)) continue;
                    if ((a.kind == kPick && b.kind == kDeliv) ||
                        (a.kind == kEmpty && b.kind == kPick) ||
                        (a.kind == kDeliv && b.kind == kEmpty)) {
                        std::swap(ops[i], ops[i + 1U]);
                        changed = true;
                    }
                }
                if (!changed) break;
            }
        }
        std::vector<std::size_t> pnext(I.pa.size(), 0U), dnext(I.dw.size(), 0U);
        std::unordered_map<int, int> traj_req;
        out.clear();
        for (std::size_t r = 0; r < rs.size(); ++r) {
            CpActionRoute route;
            route.vehicle_index = static_cast<int>(r);
            for (const Op& o : rs[r]) {
                CpActionVisit v;
                if (o.kind == kPick) {
                    auto& idx = pnext[static_cast<std::size_t>(o.cls)];
                    if (idx >= I.pmem[static_cast<std::size_t>(o.cls)].size()) { err = "pickup class overused"; return false; }
                    v.kind = CpActionKind::Pickup;
                    v.request_index = I.pmem[static_cast<std::size_t>(o.cls)][idx++];
                    traj_req[o.traj] = v.request_index;
                } else if (o.kind == kEmpty) {
                    v.kind = CpActionKind::Treatment;
                    auto it = traj_req.find(o.traj);
                    if (it == traj_req.end()) { err = "treatment before pickup"; return false; }
                    v.request_index = it->second;
                } else {
                    auto& idx = dnext[static_cast<std::size_t>(o.cls)];
                    if (idx >= I.dmem[static_cast<std::size_t>(o.cls)].size()) { err = "delivery class overused"; return false; }
                    v.kind = CpActionKind::Delivery;
                    v.request_index = I.dmem[static_cast<std::size_t>(o.cls)][idx++];
                }
                route.actions.push_back(v);
            }
            // Symmetry 40a/40c: consecutive same-location pickups (deliveries)
            // must appear in increasing request order.
            auto& acts = route.actions;
            for (std::size_t i = 0; i < acts.size();) {
                if (acts[i].kind == CpActionKind::Treatment) { ++i; continue; }
                const CpActionKind kind = acts[i].kind;
                const int loc = data.requests[static_cast<std::size_t>(acts[i].request_index)].from_id;
                std::size_t j = i + 1U;
                while (j < acts.size() && acts[j].kind == kind &&
                       data.requests[static_cast<std::size_t>(acts[j].request_index)].from_id == loc) {
                    ++j;
                }
                std::sort(acts.begin() + static_cast<std::ptrdiff_t>(i), acts.begin() + static_cast<std::ptrdiff_t>(j),
                    [](const CpActionVisit& x, const CpActionVisit& y) { return x.request_index < y.request_index; });
                i = j;
            }
            for (std::size_t i = 0; i < acts.size(); ++i) acts[i].position = static_cast<int>(i);
            out.push_back(std::move(route));
        }
        return true;
    }

    PalnsFixedKResult solve(int K, double time_limit) {
        PalnsFixedKResult res;
        t0 = Clock::now();
        has_deadline = time_limit > 0.0;
        if (has_deadline) {
            deadline = t0 + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(time_limit));
        }
        std::mt19937_64 rng(static_cast<std::uint64_t>(opt.seed) * 1000003ULL + static_cast<std::uint64_t>(K));
        next_traj = 0;
        SearchState st;
        phase = 1;

        Solution best0;
        bool have0 = false;
        const int starts = std::max(1, opt.construct_starts);
        for (int sidx = 0; sidx < starts; ++sidx) {
            if (sidx > 0 && timed_out()) break;
            Solution s = construct(K, sidx, rng);
            ++res.construction_starts_done;
            record(st, s);
            if (!have0 || s.V < best0.V - kEps ||
                (std::fabs(s.V - best0.V) <= kEps && s.C < best0.C - kEps)) {
                best0 = s;
                have0 = true;
            }
            if (best0.V <= kEps) break;
        }
        res.construction_seconds = elapsed();
        res.construction_feasible = st.have_feasible;

        if (!st.have_feasible) {
            res.phase1_iterations = run_phase(1, best0, st, rng);
        }
        // A feasible solution found by the rescue is returned without Phase II.
        const bool found_by_rescue = hub != nullptr && hub->rescue && !v2.rescue_p2;
        if (st.have_feasible && opt.phase2_enable && !found_by_rescue && !timed_out()) {
            res.phase2_iterations = run_phase(2, st.best_feasible, st, rng);
        }
        res.runtime_seconds = elapsed();
        res.time_to_feasible = st.time_to_feasible;
        res.best_violation = st.best_violation;
        if (timed_out()) res.stop_reason = "time-limit";
        else if (st.have_feasible && opt.phase2_enable) res.stop_reason = "phase2-noimp";
        else if (st.have_feasible) res.stop_reason = "feasible";
        else res.stop_reason = "time-limit";
        res.hit_time_limit = !st.have_feasible;
        if (st.have_feasible) {
            const Solution& b = st.best_feasible;
            res.best_duration = b.D;
            res.best_travel_cost = b.C;
            res.best_original_cost = original_cost(b);
            std::string err;
            if (build_action_routes(b, K, res.routes, err)) {
                res.feasible = true;
            } else {
                res.feasible = false;
                res.stop_reason = "output-error:" + err;
            }
        }
        if (!st.have_feasible && st.have_least) {
            const Solution& b = st.least;
            res.has_least_violation = true;
            res.least_violation = b.V;
            res.least_violation_duration = b.D;
            res.least_violation_travel_cost = b.C;
            std::ostringstream rt;
            for (const Route& r : b.routes) {
                res.least_violation_route_durations.push_back(r.d);
                rt << r.d << ':';
                for (std::size_t i = 0; i < r.ops.size(); ++i) {
                    const Op& o = r.ops[i];
                    rt << (i ? "," : "") << (o.kind == kPick ? 'P' : (o.kind == kEmpty ? 'E' : 'D'))
                       << o.cls << '@' << I.loc(o);
                }
                rt << ';';
            }
            res.least_violation_routes_text = rt.str();
        }
        std::ostringstream log;
        log << "[palns] k=" << K
            << " feasible=" << (res.feasible ? 1 : 0)
            << " stop_reason=" << res.stop_reason
            << " runtime_seconds=" << res.runtime_seconds
            << " construction_seconds=" << res.construction_seconds
            << " construction_starts=" << res.construction_starts_done
            << " construction_feasible=" << (res.construction_feasible ? 1 : 0)
            << " time_to_feasible=" << res.time_to_feasible
            << " best_violation=" << (std::isfinite(res.best_violation) ? res.best_violation : -1.0)
            << " best_original_cost=" << res.best_original_cost
            << " best_duration=" << res.best_duration
            << " phase1_iterations=" << res.phase1_iterations
            << " phase2_iterations=" << res.phase2_iterations;
        log << "\n[palns-trace] k=" << K << " improvements=" << st.trace.size() << " points=";
        const std::size_t tn = st.trace.size();
        for (std::size_t i = 0; i < tn; ++i) {
            if (tn > 40U && i % (tn / 40U + 1U) != 0U && i + 1U != tn) continue;
            log << st.trace[i].first << ':' << st.trace[i].second << ';';
        }
        res.log = log.str();
        return res;
    }
};

PalnsV2IncumbentSolver::PalnsV2IncumbentSolver(const SPDPData& data, const PalnsV2Options& options)
    : impl_(std::make_unique<Impl>(data, options.base)) {
    impl_->v2 = options;
    if (impl_->v2.recomb_tl > 0.0) {
        impl_->v2.pool_delta = std::max(impl_->v2.pool_delta, impl_->v2.recomb_delta);
    }
}

const std::vector<PoolBlock>& PalnsV2IncumbentSolver::pool_blocks() const {
    return impl_->last_pool;
}

PalnsV2IncumbentSolver::~PalnsV2IncumbentSolver() = default;

bool parse_palns_v2_config(const std::string& text, PalnsV2Options& options, std::string& error) {
    options.config_text = text;
    std::size_t pos = 0;
    while (pos < text.size()) {
        std::size_t end = text.find(',', pos);
        if (end == std::string::npos) end = text.size();
        const std::string item = text.substr(pos, end - pos);
        pos = end + 1U;
        if (item.empty()) continue;
        const std::size_t eq = item.find('=');
        if (eq == std::string::npos || eq == 0U) {
            error = "expected key=value, got '" + item + "'";
            return false;
        }
        const std::string key = item.substr(0, eq);
        const std::string value = item.substr(eq + 1U);
        try {
            if (!set_v2_option(options, key, value)) {
                error = "unknown key '" + key + "'";
                return false;
            }
        } catch (const std::exception&) {
            error = "bad value for '" + key + "': '" + value + "'";
            return false;
        }
    }
    return true;
}

double PalnsV2IncumbentSolver::pattern_generation_seconds() const {
    return impl_->pattern_seconds;
}

std::string PalnsV2IncumbentSolver::pattern_summary() const {
    std::size_t singles = 0;
    for (const auto& list : impl_->pool.singles) singles += list.size();
    std::ostringstream out;
    out << "[palns-pattern] pickup_classes=" << impl_->I.pa.size()
        << " delivery_classes=" << impl_->I.dw.size()
        << " singles=" << singles
        << " composite_raw=" << impl_->pool.raw
        << " composite_duplicates=" << impl_->pool.duplicates
        << " composite_dominated=" << impl_->pool.dominated
        << " composites=" << impl_->pool.comps.size()
        << " generation_seconds=" << impl_->pattern_seconds
        << " sim_alpha_raw_range=" << impl_->pool.sim_raw_range[0] << ',' << impl_->pool.sim_raw_range[1]
        << " sim_beta_raw_range=" << impl_->pool.sim_raw_range[2] << ',' << impl_->pool.sim_raw_range[3]
        << " kappa=" << impl_->opt.composite_kappa
        << " skeletons=" << palns_skeleton_set_name(impl_->opt.composite_skeletons);
    return out.str();
}

PalnsFixedKResult PalnsV2IncumbentSolver::solve_fixed_k(int vehicle_count, double time_limit_seconds) {
    if (vehicle_count <= 0) {
        throw std::runtime_error("PALNS requires a positive vehicle count.");
    }
    const PalnsOptions& O = impl_->opt;
    const int nthreads = std::max(1, O.threads);
    const long long E = std::max(1, O.sync_iters);
    const auto wall_start = Clock::now();
    auto since_start = [&]() { return std::chrono::duration<double>(Clock::now() - wall_start).count(); };
    impl_->last_pool.clear();

    SyncHub hub;
    hub.participants = nthreads;
    std::vector<std::unique_ptr<Impl>> workers;
    for (int i = 0; i < nthreads; ++i) {
        auto w = std::make_unique<Impl>(*impl_);
        w->opt.seed = O.seed + 7919U * static_cast<unsigned int>(i);
        w->hub = &hub;
        w->adopt_shared = (i % 2) == 0;
        if (O.portfolio_diverse) {
            static const double start_mul[3] = {1.0, 0.5, 2.0};
            w->opt.sa_start_ratio = O.sa_start_ratio * start_mul[i % 3];
            if ((i / 4) % 2 == 1) w->opt.kmax_ratio = O.kmax_ratio * 2.0;
        }
        if (O.portfolio_mixed_cycles && ((i / 2) % 2) == 0) {
            w->opt.sa_cycle_iters_phase2 = 0;
        }
        // Version 2: portfolio over annealing temperatures and selection rules.
        if (impl_->v2.temp_spread > 1.0) {
            const double f = impl_->v2.temp_spread;
            const double mul[3] = {1.0, f, 1.0 / f};
            w->opt.sa_start_ratio = O.sa_start_ratio * mul[i % 3];
        }
        if (impl_->v2.sel_cost_threads > 0) {
            w->v2.sel_cost = (i % impl_->v2.sel_cost_threads) == impl_->v2.sel_cost_threads - 1;
        }
        workers.push_back(std::move(w));
    }

    // Barrier action: merge the thread bests (min cost, then earliest
    // iteration, then smallest thread index) and evaluate the stop rule.
    long long epoch = 0;
    long long I_last = 0;
    long long I_first = -1;  // iteration of the first feasible shared best
    long long I_rescue = -1; // iteration at which the rescue started
    bool time_limit_hit = false;
    std::string stop_rule = "none";
    std::vector<std::tuple<long long, double, double>> trace;  // (I_last, seconds, original cost)
    std::vector<std::pair<long long, double>> clock;           // (I, seconds)
    hub.on_barrier = [&]() {
        ++epoch;
        const long long I = epoch * E;
        const double now = since_start();
        int pick = -1;
        for (int i = 0; i < nthreads; ++i) {
            const Impl& w = *workers[static_cast<std::size_t>(i)];
            if (!w.cand_has) continue;
            if (pick < 0) { pick = i; continue; }
            const Impl& p = *workers[static_cast<std::size_t>(pick)];
            if (w.cand_cost < p.cand_cost - kEps ||
                (std::fabs(w.cand_cost - p.cand_cost) <= kEps && w.cand_gi < p.cand_gi)) {
                pick = i;
            }
        }
        if (pick >= 0) {
            const Impl& p = *workers[static_cast<std::size_t>(pick)];
            if (!hub.has || p.cand_cost < hub.cost - kEps) {
                if (!hub.has) { I_first = p.cand_gi; hub.first_i = I; }
                hub.has = true;
                hub.cost = p.cand_cost;
                hub.solution = p.cand;
                I_last = p.cand_gi;
                trace.emplace_back(I_last, now, workers[0]->original_cost(p.cand));
            }
        }
        if (epoch % 8 == 0) clock.emplace_back(I, now);
        if (hub.rescue) {
            for (int i = 0; i < nthreads; ++i) {
                const Impl& w = *workers[static_cast<std::size_t>(i)];
                if (!w.rbest_has) continue;
                if (!hub.rhas || w.rbest.V < hub.rsol.V - kEps ||
                    (std::fabs(w.rbest.V - hub.rsol.V) <= kEps && w.rbest.D < hub.rsol.D - kEps)) {
                    hub.rsol = w.rbest;
                    hub.rhas = true;
                }
            }
        }
        if (hub.has) {
            const long long D = I - I_last;
            if (hub.rescue && !impl_->v2.rescue_p2) {
                hub.stop = true; stop_rule = "rescue-feasible";
            } else if (!O.phase2_enable) {
                hub.stop = true; stop_rule = "feasible";
            } else if (O.phase2_max_noimp > 0 && D >= O.phase2_max_noimp) {
                hub.stop = true; stop_rule = "max-noimp";
            } else if (O.stag_ratio > 0.0 && D >= O.stag_min_iters &&
                       static_cast<double>(D) >= O.stag_ratio * static_cast<double>(I)) {
                hub.stop = true; stop_rule = "stag-ratio";
            } else if (O.phase2_max_iters > 0 && I - I_first >= O.phase2_max_iters) {
                hub.stop = true; stop_rule = "phase2-limit";
            }
        } else if (O.phase1_max_iters > 0 && I >= O.phase1_max_iters) {
            if (O.rescue_max_iters <= 0) {
                hub.stop = true; stop_rule = "phase1-limit";
            } else if (!hub.rescue) {
                hub.rescue = true; I_rescue = I;
            } else if (I - I_rescue >= O.rescue_max_iters) {
                hub.stop = true; stop_rule = "rescue-limit";
            }
        }
        if (!hub.stop && time_limit_seconds > 0.0 && now >= time_limit_seconds) {
            hub.stop = true; stop_rule = "time-limit"; time_limit_hit = true;
        }
        if (hub.stop) clock.emplace_back(I, now);
    };

    std::vector<PalnsFixedKResult> results(static_cast<std::size_t>(nthreads));
    std::vector<std::thread> pool_threads;
    for (int i = 0; i < nthreads; ++i) {
        pool_threads.emplace_back([&, i]() {
            // The time limit is checked only at the barriers.
            results[static_cast<std::size_t>(i)] =
                workers[static_cast<std::size_t>(i)]->solve(vehicle_count, 0.0);
            hub_leave(hub);
        });
    }
    for (std::thread& t : pool_threads) t.join();

    int best = -1;
    for (int i = 0; i < nthreads; ++i) {
        const PalnsFixedKResult& r = results[static_cast<std::size_t>(i)];
        if (!r.feasible) continue;
        if (best < 0 || r.best_original_cost < results[static_cast<std::size_t>(best)].best_original_cost - 1e-9) best = i;
    }
    PalnsFixedKResult out = results[static_cast<std::size_t>(best >= 0 ? best : 0)];
    out.runtime_seconds = since_start();
    double ttf = -1.0, wait = 0.0;
    long long it1 = 0, it2 = 0;
    for (int i = 0; i < nthreads; ++i) {
        const PalnsFixedKResult& r = results[static_cast<std::size_t>(i)];
        if (r.time_to_feasible >= 0.0 && (ttf < 0.0 || r.time_to_feasible < ttf)) ttf = r.time_to_feasible;
        it1 += r.phase1_iterations;
        it2 += r.phase2_iterations;
        wait += workers[static_cast<std::size_t>(i)]->wait_seconds;
    }
    out.time_to_feasible = ttf;
    out.phase1_iterations = it1;
    out.phase2_iterations = it2;
    out.stop_reason = (best >= 0) ? stop_rule : (stop_rule == "none" ? "no-feasible" : stop_rule);
    out.hit_time_limit = best < 0;
    if (best < 0) out.feasible = false;
    else if (results[static_cast<std::size_t>(best)].stop_reason.rfind("output-error", 0) == 0) {
        out.stop_reason = results[static_cast<std::size_t>(best)].stop_reason;
    }
    std::ostringstream log;
    log << "[palns] k=" << vehicle_count
        << " feasible=" << (out.feasible ? 1 : 0)
        << " stop_reason=" << out.stop_reason
        << " runtime_seconds=" << out.runtime_seconds
        << " construction_seconds=" << out.construction_seconds
        << " construction_starts=" << out.construction_starts_done
        << " construction_feasible=" << (out.construction_feasible ? 1 : 0)
        << " time_to_feasible=" << out.time_to_feasible
        << " best_violation=" << (std::isfinite(out.best_violation) ? out.best_violation : -1.0)
        << " best_original_cost=" << out.best_original_cost
        << " best_duration=" << out.best_duration
        << " phase1_iterations=" << out.phase1_iterations
        << " phase2_iterations=" << out.phase2_iterations
        << " threads=" << nthreads
        << " best_thread=" << best
        << " sync_iters=" << E
        << " epochs=" << epoch
        << " iterations_per_thread=" << epoch * E
        << " last_improvement_iteration=" << I_last
        << " barrier_wait_seconds=" << wait
        << " time_limit_hit=" << (time_limit_hit ? 1 : 0)
        << " phase1_max_iters=" << O.phase1_max_iters
        << " phase2_max_iters=" << O.phase2_max_iters
        << " first_feasible_iteration=" << I_first
        << " rescue_max_iters=" << O.rescue_max_iters
        << " rescue_routes=" << O.rescue_routes
        << " rescue_start_iteration=" << I_rescue << '\n';
    if (best < 0) {
        int least = -1;
        for (int i = 0; i < nthreads; ++i) {
            const PalnsFixedKResult& r = results[static_cast<std::size_t>(i)];
            if (!r.has_least_violation) continue;
            if (least < 0) { least = i; continue; }
            const PalnsFixedKResult& l = results[static_cast<std::size_t>(least)];
            if (r.least_violation < l.least_violation - 1e-9 ||
                (std::fabs(r.least_violation - l.least_violation) <= 1e-9 &&
                 r.least_violation_travel_cost < l.least_violation_travel_cost - 1e-9)) {
                least = i;
            }
        }
        for (int i = 0; i < nthreads; ++i) {
            const PalnsFixedKResult& r = results[static_cast<std::size_t>(i)];
            if (!r.has_least_violation) continue;
            log << "[palns-least-violation] k=" << vehicle_count << " thread=" << i
                << " violation=" << r.least_violation << " duration=" << r.least_violation_duration
                << " travel_cost=" << r.least_violation_travel_cost << " route_durations=";
            for (double d : r.least_violation_route_durations) log << d << ',';
            log << '\n';
        }
        if (least >= 0) {
            log << "[palns-least-violation-routes] k=" << vehicle_count << " thread=" << least << " routes="
                << results[static_cast<std::size_t>(least)].least_violation_routes_text << '\n';
        }
    }
    for (int i = 0; i < nthreads; ++i) {
        const PalnsFixedKResult& r = results[static_cast<std::size_t>(i)];
        log << "[palns-thread] index=" << i << " feasible=" << (r.feasible ? 1 : 0)
            << " cost=" << r.best_original_cost
            << " time_to_feasible=" << r.time_to_feasible
            << " iterations=" << workers[static_cast<std::size_t>(i)]->gi
            << " rescue_start_violation=" << workers[static_cast<std::size_t>(i)]->rescue_start_V
            << " barrier_wait_seconds=" << workers[static_cast<std::size_t>(i)]->wait_seconds << '\n';
        log << "[palns-v2-thread] index=" << i
            << " focus_iters=" << workers[static_cast<std::size_t>(i)]->stat_focus_iters
            << " focus_episodes=" << workers[static_cast<std::size_t>(i)]->stat_focus_episodes
            << " focus_best=" << workers[static_cast<std::size_t>(i)]->stat_focus_best
            << " global_best=" << workers[static_cast<std::size_t>(i)]->stat_global_best
            << " reloc_hit=" << workers[static_cast<std::size_t>(i)]->stat_reloc_hit
            << " reloc_miss=" << workers[static_cast<std::size_t>(i)]->stat_reloc_miss << '\n';
        if (impl_->v2.profile) {
            const auto& pf = workers[static_cast<std::size_t>(i)]->prof;
            log << "[palns-v2-profile] index=" << i
                << " ruin_s=" << pf.ruin << " ruin_n=" << pf.n_ruin
                << " dswap_s=" << pf.dswap << " dswap_n=" << pf.n_dswap
                << " blocks_s=" << pf.blocks << " blocks_n=" << pf.n_blocks
                << " reloc_s=" << pf.reloc << " reloc_n=" << pf.n_reloc
                << " ls_calls=" << pf.ls_calls << " ls_moves=" << pf.ls_moves
                << " wait_s=" << workers[static_cast<std::size_t>(i)]->wait_seconds << '\n';
        }
        if (impl_->v2.diag) {
            const auto& dg = workers[static_cast<std::size_t>(i)]->diag;
            auto arr = [&](const char* name, const long long* v, int n) {
                log << ' ' << name << '=';
                for (int x = 0; x < n; ++x) log << (x ? "," : "") << v[x];
            };
            log << "[palns-v2-diag] index=" << i
                << " iters=" << dg.iters << " feasible=" << dg.feasible << " seen_before=" << dg.seen_before
                << " accepted=" << dg.accepted << " uphill=" << dg.uphill
                << " comp_iters=" << dg.comp_iters << " comp_cand=" << dg.comp_cand << " comp_ins=" << dg.comp_ins
                << " single_ins=" << dg.single_ins << " single_split=" << dg.single_split
                << " ls_dswap=" << dg.ls_dswap << " ls_block_swap=" << dg.ls_blocks[0]
                << " ls_block_relocate=" << dg.ls_blocks[1] << " ls_block_move=" << dg.ls_blocks[2]
                << " ls_traj_relocate=" << dg.ls_reloc
                << " best_after_repair=" << dg.best_after_repair << " best_after_light=" << dg.best_after_light
                << " best_after_deep=" << dg.best_after_deep << " lam_back=" << dg.lam_back
                << " indep_adopt=" << dg.indep_adopt;
            arr("d_use", dg.d_use, 5); arr("d_acc", dg.d_acc, 5); arr("d_imp", dg.d_imp, 5); arr("d_best", dg.d_best, 5);
            arr("q_use", dg.q_use, 3); arr("q_acc", dg.q_acc, 3); arr("q_imp", dg.q_imp, 3); arr("q_best", dg.q_best, 3);
            arr("r_use", dg.r_use, 5); arr("r_acc", dg.r_acc, 5); arr("r_imp", dg.r_imp, 5); arr("r_best", dg.r_best, 5);
            arr("k_use", dg.k_use, 16); arr("k_best", dg.k_best, 16);
            log << '\n' << "[palns-v2-diag-cycles] index=" << i << " cycles=";
            for (const auto& cy : dg.cycles) {
                log << cy.start << ':' << (std::isfinite(cy.low) ? cy.low : -1.0) << ':' << cy.n << ':' << cy.acc << ':'
                    << cy.feas << ':' << cy.best << ':' << (cy.n > 0 ? cy.gap_sum / static_cast<double>(cy.n) : 0.0)
                    << ':' << cy.lam << ':' << cy.cur_inf << ':'
                    << (cy.cur_inf > 0 ? cy.cur_inf_V / static_cast<double>(cy.cur_inf) : 0.0) << ':'
                    << (std::isfinite(cy.lam_lo) ? cy.lam_lo : 0.0) << ':' << cy.lam_hi << ':' << cy.mul << ';';
            }
            log << '\n';
        }
    }
    log << "[palns-v2] config=" << (impl_->v2.config_text.empty() ? "(none)" : impl_->v2.config_text) << '\n';
    if (impl_->v2.pool_ratio > 0.0 || impl_->v2.pool_delta > 0.0) {
        // Pool for the recombination step: the blocks of all threads with the
        // cheapest solution that contained each, as (kind, location, type).
        {
            const Instance& I = impl_->I;
            std::unordered_map<std::string, double> merged;
            for (int i = 0; i < nthreads; ++i) {
                const Impl& w = *workers[static_cast<std::size_t>(i)];
                for (std::size_t b = 0; b < w.pool_keys.size(); ++b) {
                    auto it = merged.find(w.pool_keys[b]);
                    if (it == merged.end()) merged.emplace(w.pool_keys[b], w.pool_min[b]);
                    else if (w.pool_min[b] < it->second) it->second = w.pool_min[b];
                }
            }
            std::vector<const std::pair<const std::string, double>*> order;
            order.reserve(merged.size());
            for (const auto& entry : merged) order.push_back(&entry);
            std::sort(order.begin(), order.end(),
                      [](const auto* a, const auto* b) { return a->first < b->first; });
            impl_->last_pool.reserve(order.size());
            for (const auto* entry : order) {
                const std::string& key = entry->first;
                PoolBlock block;
                block.min_cost = entry->second;
                block.ops.reserve(key.size() / 3U);
                for (std::size_t x = 0; x + 2U < key.size(); x += 3U) {
                    const int kind = static_cast<unsigned char>(key[x]);
                    const std::size_t cls = static_cast<std::size_t>(
                        static_cast<unsigned char>(key[x + 1U]) |
                        (static_cast<unsigned char>(key[x + 2U]) << 8));
                    // Container type as in the instance data (I.pf / I.df are compact indices).
                    const int member = kind == kDeliv ? I.dmem[cls].front() : I.pmem[cls].front();
                    const int type = impl_->data.requests[static_cast<std::size_t>(member)].container_type;
                    if (kind == kPick) block.ops.push_back({'P', I.pa[cls], type});
                    else if (kind == kEmpty) block.ops.push_back({'E', I.pb[cls], type});
                    else block.ops.push_back({'D', I.dw[cls], type});
                }
                impl_->last_pool.push_back(std::move(block));
            }
            log << "[palns-v2-pool] ratio=" << impl_->v2.pool_ratio << " delta=" << impl_->v2.pool_delta
                << " merged_blocks=" << impl_->last_pool.size() << '\n';
        }
        // Block pool (diagnostics). Lines of the file, tab separated:
        //   K  vehicles requests threads ratio
        //   PC class location facility type requests...   (pickup class; requests 0-based)
        //   DC class location type requests...            (delivery class)
        //   TR iteration seconds cost                     (shared best)
        //   CL iteration seconds                          (iteration clock)
        //   B  thread block ops...                        (P<class>, E<class>, D<class>)
        //   V  thread block iteration cost                (cheapest solution containing the block so far)
        const Instance& I = impl_->I;
        std::ofstream pf;
        if (!impl_->v2.pool_file.empty()) pf.open(impl_->v2.pool_file, std::ios::app);
        if (pf) {
            pf.precision(10);
            pf << "K\t" << vehicle_count << '\t' << I.n << '\t' << nthreads << '\t' << impl_->v2.pool_ratio << '\n';
            for (std::size_t p = 0; p < I.pa.size(); ++p) {
                pf << "PC\t" << p << '\t' << I.pa[p] << '\t' << I.pb[p] << '\t' << I.pf[p];
                for (int r : I.pmem[p]) pf << '\t' << r;
                pf << '\n';
            }
            for (std::size_t q = 0; q < I.dw.size(); ++q) {
                pf << "DC\t" << q << '\t' << I.dw[q] << '\t' << I.df[q];
                for (int r : I.dmem[q]) pf << '\t' << r;
                pf << '\n';
            }
            for (const auto& p : trace) pf << "TR\t" << std::get<0>(p) << '\t' << std::get<1>(p) << '\t' << std::get<2>(p) << '\n';
            for (const auto& p : clock) pf << "CL\t" << p.first << '\t' << p.second << '\n';
        }
        for (int i = 0; i < nthreads; ++i) {
            const Impl& w = *workers[static_cast<std::size_t>(i)];
            log << "[palns-v2-pool] index=" << i << " solutions=" << w.pool_solutions
                << " blocks=" << w.pool_keys.size() << " events=" << w.pool_events.size() << '\n';
            if (!pf) continue;
            for (std::size_t b = 0; b < w.pool_keys.size(); ++b) {
                const std::string& key = w.pool_keys[b];
                pf << "B\t" << i << '\t' << b;
                for (std::size_t x = 0; x + 2U < key.size(); x += 3U) {
                    const int kind = static_cast<unsigned char>(key[x]);
                    const int cls = static_cast<unsigned char>(key[x + 1U]) |
                        (static_cast<unsigned char>(key[x + 2U]) << 8);
                    pf << '\t' << (kind == kPick ? 'P' : (kind == kEmpty ? 'E' : 'D')) << cls;
                }
                pf << '\n';
            }
            for (const auto& e : w.pool_events) {
                pf << "V\t" << i << '\t' << e.block << '\t' << e.gi << '\t' << e.cost << '\n';
            }
        }
        if (!impl_->v2.pool_file.empty() && !pf) log << "[palns-v2-pool] cannot write " << impl_->v2.pool_file << '\n';
    }
    log << "[palns-shared-trace] k=" << vehicle_count << " improvements=" << trace.size() << " points=";
    for (const auto& p : trace) log << std::get<0>(p) << ':' << std::get<1>(p) << ':' << std::get<2>(p) << ';';
    log << "\n[palns-iter-clock] k=" << vehicle_count << " points=";
    for (const auto& p : clock) log << p.first << ':' << p.second << ';';
    out.log = log.str();
    return out;
}

}  // namespace spdp

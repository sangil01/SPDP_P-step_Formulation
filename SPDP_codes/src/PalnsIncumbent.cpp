#include "PalnsIncumbent.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <thread>
#include <chrono>
#include <cmath>
#include <cstdint>
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

const char* palns_skeleton_set_name(PalnsSkeletonSet value) {
    return value == PalnsSkeletonSet::All ? "all" : "S3,S5";
}
const char* palns_composite_filter_name(PalnsCompositeFilter value) {
    return value == PalnsCompositeFilter::DeltaPositive ? "delta-positive" : "all";
}
const char* palns_construct_order_name(PalnsConstructOrder value) {
    switch (value) {
    case PalnsConstructOrder::Lpt: return "lpt";
    case PalnsConstructOrder::Random: return "random";
    case PalnsConstructOrder::Regret: return "regret";
    case PalnsConstructOrder::Mixed: return "mixed";
    }
    return "mixed";
}
const char* palns_ls_frequency_name(PalnsLsFrequency value) {
    switch (value) {
    case PalnsLsFrequency::Every: return "every";
    case PalnsLsFrequency::OnImprove: return "on-improve";
    case PalnsLsFrequency::EveryM: return "every-m";
    }
    return "every";
}
bool parse_palns_skeleton_set(const std::string& text, PalnsSkeletonSet& value) {
    if (text == "all") { value = PalnsSkeletonSet::All; return true; }
    if (text == "S3,S5" || text == "s3,s5" || text == "S3S5") {
        value = PalnsSkeletonSet::S3S5; return true;
    }
    return false;
}
bool parse_palns_composite_filter(const std::string& text, PalnsCompositeFilter& value) {
    if (text == "delta-positive") { value = PalnsCompositeFilter::DeltaPositive; return true; }
    if (text == "all") { value = PalnsCompositeFilter::All; return true; }
    return false;
}
bool parse_palns_construct_order(const std::string& text, PalnsConstructOrder& value) {
    if (text == "lpt") { value = PalnsConstructOrder::Lpt; return true; }
    if (text == "random") { value = PalnsConstructOrder::Random; return true; }
    if (text == "regret") { value = PalnsConstructOrder::Regret; return true; }
    if (text == "mixed") { value = PalnsConstructOrder::Mixed; return true; }
    return false;
}
bool parse_palns_ls_frequency(const std::string& text, PalnsLsFrequency& value) {
    if (text == "every") { value = PalnsLsFrequency::Every; return true; }
    if (text == "on-improve") { value = PalnsLsFrequency::OnImprove; return true; }
    if (text == "every-m") { value = PalnsLsFrequency::EveryM; return true; }
    return false;
}

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

    pool.sim.assign(static_cast<std::size_t>(P), std::vector<double>(static_cast<std::size_t>(P), 0.0));
    for (int p = 0; p < P; ++p) {
        for (int r = 0; r < P; ++r) {
            const double s1 = 1.0 / (1.0 + I.tt(I.pa[static_cast<std::size_t>(p)], I.pa[static_cast<std::size_t>(r)]));
            const double s2 = 1.0 / (1.0 + I.tt(I.pb[static_cast<std::size_t>(p)], I.pb[static_cast<std::size_t>(r)]));
            const double s3 = I.pf[static_cast<std::size_t>(p)] == I.pf[static_cast<std::size_t>(r)] ? 1.0 : 0.0;
            pool.sim[static_cast<std::size_t>(p)][static_cast<std::size_t>(r)] =
                opt.sim_s1 * s1 + opt.sim_s2 * s2 + opt.sim_s3 * s3;
        }
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
std::vector<Trip> assign_singles(const Instance& I, const Counts& k, int mode, std::mt19937_64& rng) {
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
                    cost[i][j] = I.single_d(U[i], Z[j]) + 1e-6 * I.single_c(U[i], Z[j]);
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
                        w[j] = 1.0 / std::max(1e-9, I.single_d(p, rem[j]));
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

// Evaluates the best insertion of block B into every route. per_route (if
// given) receives the best position per route (for regret).
InsEval eval_insertion(const Instance& I, const Solution& s, const BlockInfo& B, bool allow_split,
                       bool phase2, double lambda, std::vector<InsEval>* per_route) {
    InsEval best;
    if (per_route != nullptr) per_route->assign(s.routes.size(), InsEval{});
    RouteCtx ctx;
    for (std::size_t r = 0; r < s.routes.size(); ++r) {
        const Route& route = s.routes[r];
        make_ctx(I, route, ctx);
        const int m = static_cast<int>(route.ops.size());
        InsEval rbest;
        auto offer = [&](InsEval e) {
            e.route = static_cast<int>(r);
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
            for (int pass = 0; pass < (phase2 ? 2 : 1); ++pass) {
                int bi, bj, bk;
                double bv;
                best_split_by(I, ctx, B, pass == 1, bi, bj, bk, bv);
                if (bi < 0) continue;
                InsEval e;
                e.split = true;
                e.i = bi; e.j = bj; e.k = bk;
                exact_split(I, ctx, B, bi, bj, bk, e.dd, e.dc);
                offer(e);
            }
        }
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

}  // namespace

// ---------------------------------------------------------------- solver

namespace {
// Best feasible solution shared by the portfolio threads.
struct SharedBest {
    std::mutex mutex;
    bool has = false;
    double cost = 0.0;
    Solution solution;
};
}  // namespace

struct PalnsIncumbentSolver::Impl {
    const SPDPData& data;
    PalnsOptions opt;
    Instance I;
    PatternPool pool;
    double pattern_seconds = 0.0;
    int next_traj = 0;
    int phase = 1;
    double lambda = 1.0;  // Phase II penalty on duration excess
    Clock::time_point t0;
    Clock::time_point deadline;
    bool has_deadline = false;
    std::atomic<bool>* stop_flag = nullptr;
    SharedBest* shared = nullptr;

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

    void publish(const Solution& s) {
        if (shared == nullptr || s.V > kEps) return;
        std::lock_guard<std::mutex> lock(shared->mutex);
        if (!shared->has || s.C < shared->cost - kEps) {
            shared->has = true;
            shared->cost = s.C;
            shared->solution = s;
        }
    }

    bool adopt(Solution& out, double own_cost) {
        if (shared == nullptr) return false;
        std::lock_guard<std::mutex> lock(shared->mutex);
        if (!shared->has || shared->cost >= own_cost - kEps) return false;
        out = shared->solution;
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
        if (stop_flag != nullptr && stop_flag->load(std::memory_order_relaxed)) return true;
        return has_deadline && Clock::now() >= deadline;
    }

    double objective(const Solution& s) const {
        return phase == 1 ? s.V + opt.eta * s.D : s.C + lambda * s.V;
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
        if (mode == 0) {
            int best = ok.front();
            for (int idx : ok) {
                const Composite& a = pool.comps[static_cast<std::size_t>(idx)];
                const Composite& b = pool.comps[static_cast<std::size_t>(best)];
                if (a.d < b.d - kEps || (std::fabs(a.d - b.d) <= kEps && a.c < b.c - kEps)) best = idx;
            }
            return best;
        }
        if (mode == 1) {
            std::uniform_int_distribution<std::size_t> dist(0, ok.size() - 1U);
            return ok[dist(rng)];
        }
        std::vector<double> w(ok.size());
        for (std::size_t i = 0; i < ok.size(); ++i) {
            w[i] = 1.0 / std::max(1e-9, pool.comps[static_cast<std::size_t>(ok[i])].d);
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
        std::vector<InsEval> per_route;
        while (true) {
            int pick = -1;
            InsEval pick_eval;
            double pk1 = kInf, pk2 = kInf, pk3 = kInf;  // ordering keys
            for (std::size_t x = 0; x < cand.size(); ++x) {
                if (done[x]) continue;
                if (cnt != nullptr && !trip_valid(cand[x], *cnt)) { done[x] = 1; continue; }
                const InsEval e = eval_insertion(I, s, blocks[x], allow_split, phase2, lambda,
                                                 regret ? &per_route : nullptr);
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
            apply_insertion(I, s, blocks[static_cast<std::size_t>(pick)], pick_eval);
            if (cnt != nullptr) consume_trip(cand[static_cast<std::size_t>(pick)], *cnt);
            done[static_cast<std::size_t>(pick)] = 1;
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
            return (excess(r.d, I.T) - excess(d, I.T)) + opt.eta * (r.d - d);
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
        std::vector<Trip> singles = assign_singles(I, cnt, qop, mrng);
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
            std::vector<Trip> rest = assign_singles(I, cnt, 0, rng);
            std::stable_sort(rest.begin(), rest.end(), [](const Trip& a, const Trip& b) { return a.d > b.d; });
            for (const Trip& t : rest) insert_trip(s, t, true, phase2);
        }
        refresh_totals(I, s);
    }

    // ------------------------------------------------------------ local search
    // Objective change of replacing route values (d,c) by (nd,nc).
    double route_delta(double d, double c, double nd, double nc) const {
        if (phase == 1) {
            return (excess(nd, I.T) - excess(d, I.T)) + opt.eta * (nd - d);
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
                f = phase == 1 ? V + opt.eta * D : C + lambda * V;
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

    void local_search(Solution& s, bool deep) {
        int moves = 0;
        bool improved = true;
        while (improved && moves < 200 && !timed_out()) {
            improved = false;
            if (ls_dswap(s) || ls_blocks(s, 0) || ls_blocks(s, 1) || ls_blocks(s, 2) ||
                (deep && ls_relocate(s))) {
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
        double time_to_target = -1.0;
        bool reached_target = false;
        std::vector<std::pair<double, double>> trace;  // (time, original cost)
    };

    void record(SearchState& st, const Solution& s) {
        st.best_violation = std::min(st.best_violation, s.V);
        if (s.V > kEps) return;
        if (!st.have_feasible || s.C < st.best_feasible.C - kEps) {
            st.best_feasible = s;
            st.have_feasible = true;
            st.trace.push_back({elapsed(), original_cost(s)});
        }
        if (st.time_to_feasible < 0.0) st.time_to_feasible = elapsed();
        if (opt.target_cost >= 0.0 && !st.reached_target &&
            original_cost(st.best_feasible) <= opt.target_cost + 1e-6) {
            st.reached_target = true;
            st.time_to_target = elapsed();
        }
    }

    long long run_phase(int which, Solution start, SearchState& st, std::mt19937_64& rng) {
        phase = which;
        refresh_totals(I, start);
        if (which == 2) {
            lambda = std::max(1.0, 2.0 * start.C / std::max(1.0, start.D));
        }
        Solution S = start;
        Solution best = start;
        double fS = objective(S);
        double fBest = which == 1 ? fS : start.C;
        int lam_window = 0, lam_feasible = 0;
        std::vector<double> wd(5, 1.0), wq(3, 1.0), wr(5, 1.0);
        std::vector<double> sd(5, 0.0), sq(3, 0.0), sr(5, 0.0);
        std::vector<int> nd(5, 0), nq(3, 0), nr(5, 0);
        const int weight_period = 50;
        const double rho = 0.1;
        const double theta0 = std::max(1e-9, opt.sa_start_ratio * std::max(fS, 1e-9) / std::log(2.0));
        const long long cycle = std::max<long long>(1, opt.sa_cycle_iters);
        std::unordered_set<std::uint64_t> seen;
        seen.insert(solution_hash(S));
        const int n = I.n;
        const int kmax = std::max(1, std::min(n, std::max(opt.kmax_min,
            static_cast<int>(std::ceil(opt.kmax_ratio * static_cast<double>(n))))));
        std::uniform_int_distribution<int> kdist(1, kmax);
        std::uniform_real_distribution<double> U01(0.0, 1.0);
        long long iter = 0;
        long long noimp = 0;
        while (!timed_out()) {
            if (which == 1 && st.have_feasible) break;
            if (st.reached_target) break;
            if (which == 2 && noimp >= opt.phase2_max_noimp) break;
            ++iter;
            std::vector<double> wdd = wd;
            if (which == 2) wdd[4] = 0.0;
            const int dop = roulette_index(wdd, rng);
            const int qop = roulette_index(wq, rng);
            const int rop = roulette_index(wr, rng);
            const int k = kdist(rng);
            Solution sp = S;
            int ref_class = -1;
            const std::unordered_set<int> rem = destroy(sp, dop, k, rng, ref_class);
            Counts cnt{std::vector<int>(static_cast<std::size_t>(npc()), 0),
                       std::vector<int>(static_cast<std::size_t>(ndc()), 0)};
            remove_trajs(I, sp, rem, cnt);
            std::vector<Trip> cand = select_candidates(cnt, qop, rng);
            repair(sp, std::move(cand), cnt, rop, ref_class, rng);
            bool run_ls = true;
            if (opt.ls_frequency == PalnsLsFrequency::OnImprove) run_ls = objective(sp) < fS - 1e-9;
            else if (opt.ls_frequency == PalnsLsFrequency::EveryM) run_ls = (iter % std::max(1, opt.ls_every_m)) == 0;
            if (run_ls) {
                local_search(sp, false);
                const double fl = objective(sp);
                const double ref = which == 1 ? fBest : fBest;
                if (fl <= ref + opt.ls_deep_ratio * std::fabs(ref) + 1e-9) local_search(sp, true);
            }
            refresh_totals(I, sp);
            record(st, sp);
            double reward = 0.0;
            const double fsp = objective(sp);
            {
                const double delta = fsp - fS;
                bool accept = false;
                if (delta < -1e-12) {
                    accept = true;
                } else {
                    const double pos = static_cast<double>((iter - 1) % cycle) / static_cast<double>(cycle);
                    const double theta = theta0 * std::pow(opt.sa_end_ratio, pos);
                    accept = U01(rng) < std::exp(-delta / std::max(theta, 1e-300));
                }
                const std::uint64_t h = solution_hash(sp);
                const bool is_new = seen.find(h) == seen.end();
                const bool new_best = which == 1 ? (fsp < fBest - 1e-9)
                                                 : (sp.V <= kEps && sp.C < fBest - 1e-9);
                if (new_best) {
                    best = sp;
                    fBest = which == 1 ? fsp : sp.C;
                    noimp = 0;
                    reward = 10.0;
                    if (which == 2) publish(sp);
                } else {
                    ++noimp;
                    if (fsp < fS - 1e-9) reward = 3.0;
                    else if (accept && is_new) reward = 1.0;
                }
                if (accept) {
                    S = sp;
                    fS = fsp;
                    seen.insert(h);
                }
            }
            if (which == 2) {
                ++lam_window;
                if (sp.V <= kEps) ++lam_feasible;
                if (lam_window >= 100) {
                    const double frac = static_cast<double>(lam_feasible) / lam_window;
                    if (frac < 0.4) lambda = std::min(1e6, lambda * 1.5);
                    else if (frac > 0.6) lambda = std::max(1e-3, lambda / 1.5);
                    lam_window = 0;
                    lam_feasible = 0;
                    fS = objective(S);
                }
            }
            if (opt.sa_restart_best && iter % cycle == 0) {
                if (which == 2) {
                    Solution g;
                    if (adopt(g, fBest)) {
                        best = g;
                        fBest = g.C;
                        record(st, g);
                    }
                }
                S = best;
                fS = objective(S);
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
        }
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

        if (!st.have_feasible && !st.reached_target) {
            res.phase1_iterations = run_phase(1, best0, st, rng);
        }
        if (st.have_feasible && !st.reached_target && opt.phase2_enable && !timed_out()) {
            res.phase2_iterations = run_phase(2, st.best_feasible, st, rng);
        }
        res.runtime_seconds = elapsed();
        res.time_to_feasible = st.time_to_feasible;
        res.time_to_target = st.time_to_target;
        res.reached_target = st.reached_target;
        res.best_violation = st.best_violation;
        if (st.reached_target) res.stop_reason = "target";
        else if (timed_out()) res.stop_reason = "time-limit";
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
        std::ostringstream log;
        log << "[palns] k=" << K
            << " feasible=" << (res.feasible ? 1 : 0)
            << " stop_reason=" << res.stop_reason
            << " runtime_seconds=" << res.runtime_seconds
            << " construction_seconds=" << res.construction_seconds
            << " construction_starts=" << res.construction_starts_done
            << " construction_feasible=" << (res.construction_feasible ? 1 : 0)
            << " time_to_feasible=" << res.time_to_feasible
            << " time_to_target=" << res.time_to_target
            << " target_cost=" << opt.target_cost
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

PalnsIncumbentSolver::PalnsIncumbentSolver(const SPDPData& data, const PalnsOptions& options)
    : impl_(std::make_unique<Impl>(data, options)) {}

PalnsIncumbentSolver::~PalnsIncumbentSolver() = default;

double PalnsIncumbentSolver::pattern_generation_seconds() const {
    return impl_->pattern_seconds;
}

std::string PalnsIncumbentSolver::pattern_summary() const {
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
        << " kappa=" << impl_->opt.composite_kappa
        << " skeletons=" << palns_skeleton_set_name(impl_->opt.composite_skeletons);
    return out.str();
}

PalnsFixedKResult PalnsIncumbentSolver::solve_fixed_k(int vehicle_count, double time_limit_seconds) {
    if (vehicle_count <= 0) {
        throw std::runtime_error("PALNS requires a positive vehicle count.");
    }
    const int nthreads = std::max(1, impl_->opt.threads);
    if (nthreads == 1) {
        return impl_->solve(vehicle_count, time_limit_seconds);
    }
    const auto wall_start = Clock::now();
    std::atomic<bool> stop{false};
    SharedBest shared;
    std::vector<std::unique_ptr<Impl>> workers;
    for (int i = 0; i < nthreads; ++i) {
        auto w = std::make_unique<Impl>(*impl_);
        w->opt.seed = impl_->opt.seed + 7919U * static_cast<unsigned int>(i);
        w->stop_flag = &stop;
        w->shared = &shared;
        workers.push_back(std::move(w));
    }
    std::vector<PalnsFixedKResult> results(static_cast<std::size_t>(nthreads));
    std::vector<std::thread> pool_threads;
    for (int i = 0; i < nthreads; ++i) {
        pool_threads.emplace_back([&, i]() {
            results[static_cast<std::size_t>(i)] =
                workers[static_cast<std::size_t>(i)]->solve(vehicle_count, time_limit_seconds);
            stop.store(true);
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
    out.runtime_seconds = std::chrono::duration<double>(Clock::now() - wall_start).count();
    double ttf = -1.0, ttt = -1.0;
    long long it1 = 0, it2 = 0;
    for (const PalnsFixedKResult& r : results) {
        if (r.time_to_feasible >= 0.0 && (ttf < 0.0 || r.time_to_feasible < ttf)) ttf = r.time_to_feasible;
        if (r.time_to_target >= 0.0 && (ttt < 0.0 || r.time_to_target < ttt)) ttt = r.time_to_target;
        it1 += r.phase1_iterations;
        it2 += r.phase2_iterations;
        out.reached_target = out.reached_target || r.reached_target;
    }
    out.time_to_feasible = ttf;
    out.time_to_target = ttt;
    out.phase1_iterations = it1;
    out.phase2_iterations = it2;
    if (best < 0) {
        out.feasible = false;
        out.hit_time_limit = true;
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
        << " time_to_target=" << out.time_to_target
        << " target_cost=" << impl_->opt.target_cost
        << " best_violation=" << (std::isfinite(out.best_violation) ? out.best_violation : -1.0)
        << " best_original_cost=" << out.best_original_cost
        << " best_duration=" << out.best_duration
        << " phase1_iterations=" << out.phase1_iterations
        << " phase2_iterations=" << out.phase2_iterations
        << " threads=" << nthreads
        << " best_thread=" << best << '\n';
    for (int i = 0; i < nthreads; ++i) {
        const PalnsFixedKResult& r = results[static_cast<std::size_t>(i)];
        log << "[palns-thread] index=" << i << " feasible=" << (r.feasible ? 1 : 0)
            << " stop_reason=" << r.stop_reason << " cost=" << r.best_original_cost
            << " time_to_feasible=" << r.time_to_feasible << " time_to_target=" << r.time_to_target << '\n';
    }
    out.log = log.str();
    if (!out.log.empty() && out.log.back() == '\n') out.log.pop_back();
    return out;
}

// ---------------------------------------------------------------- mapping

PalnsMappedIncumbent map_palns_routes_to_multigraph(
    const SPDPData& data,
    const MultiDiGraph& graph,
    const std::vector<CpActionRoute>& routes,
    int expected_vehicle_count
) {
    PalnsMappedIncumbent mapped;
    try {
        const int n = static_cast<int>(data.requests.size());
        std::set<int> seen_nodes;
        for (const CpActionRoute& route : routes) {
            RouteState state;
            NodeId current_node = 0;
            int current_location = -1;
            State start_state = state.canonical_state();
            std::vector<int> treatments;
            std::vector<int> treatment_requests;
            int last_treatment = -1;

            auto flush_to = [&](NodeId next_node, int next_location, double destination_service,
                                const CpActionVisit* destination) {
                for (int request_index : treatment_requests) {
                    const Request& tr = data.requests[static_cast<std::size_t>(request_index)];
                    if (!state.try_empty_request(request_index, tr.to_id)) {
                        throw std::runtime_error("PALNS route empties a skip that is not onboard.");
                    }
                }
                if (destination != nullptr) {
                    const Request& req = data.requests[static_cast<std::size_t>(destination->request_index)];
                    if (destination->kind == CpActionKind::Pickup) {
                        if (!state.try_pickup(OnboardSkip{destination->request_index, req.container_type, req.to_id, true})) {
                            throw std::runtime_error("PALNS route exceeds capacity.");
                        }
                    } else if (!state.try_delivery(req.container_type)) {
                        throw std::runtime_error("PALNS route delivers without a compatible empty skip.");
                    }
                }
                const State end_state = state.canonical_state();
                double travel = 0.0, cost = 0.0;
                int cur = current_location;
                for (int l : treatments) {
                    if (cur >= 0) {
                        travel += data.time[static_cast<std::size_t>(cur)][static_cast<std::size_t>(l)];
                        cost += data.distance[static_cast<std::size_t>(cur)][static_cast<std::size_t>(l)];
                    }
                    cur = l;
                }
                if (cur >= 0 && next_location >= 0) {
                    travel += data.time[static_cast<std::size_t>(cur)][static_cast<std::size_t>(next_location)];
                    cost += data.distance[static_cast<std::size_t>(cur)][static_cast<std::size_t>(next_location)];
                }
                const double our_time = travel + static_cast<double>(treatment_requests.size()) * data.time_empty + destination_service;
                if (current_node == 0 && destination != nullptr) cost += data.fixed_vehicle_cost;
                int exact = -1, dominating = -1, fastest = -1;
                for (std::size_t edge_id : graph.outgoing_edge_indices(current_node)) {
                    const EdgeRecord& e = graph.edges()[edge_id];
                    if (e.v != next_node || e.data.start_state != start_state || e.data.end_state != end_state) continue;
                    const int id = static_cast<int>(edge_id);
                    if (e.data.sequence_pi == treatments && std::fabs(e.data.time - our_time) <= 1e-6) {
                        if (exact < 0 || e.data.cost < graph.edges()[static_cast<std::size_t>(exact)].data.cost) exact = id;
                    }
                    if (e.data.time <= our_time + 1e-6 && e.data.cost <= cost + 1e-6) {
                        if (dominating < 0 || e.data.time < graph.edges()[static_cast<std::size_t>(dominating)].data.time) dominating = id;
                    }
                    if (fastest < 0 || e.data.time < graph.edges()[static_cast<std::size_t>(fastest)].data.time ||
                        (e.data.time == graph.edges()[static_cast<std::size_t>(fastest)].data.time &&
                         e.data.cost < graph.edges()[static_cast<std::size_t>(fastest)].data.cost)) {
                        fastest = id;
                    }
                }
                int chosen = exact >= 0 ? exact : (dominating >= 0 ? dominating : fastest);
                if (chosen < 0) {
                    std::ostringstream error;
                    error << "No multigraph edge for PALNS transition " << current_node << " -> " << next_node << '.';
                    throw std::runtime_error(error.str());
                }
                mapped.active_edge_ids.push_back(chosen);
                mapped.total_duration += graph.edges()[static_cast<std::size_t>(chosen)].data.time;
                mapped.total_original_cost += graph.edges()[static_cast<std::size_t>(chosen)].data.cost;
                current_node = next_node;
                current_location = next_location;
                start_state = end_state;
                treatments.clear();
                treatment_requests.clear();
            };

            for (const CpActionVisit& visit : route.actions) {
                if (visit.request_index < 0 || visit.request_index >= n) {
                    throw std::runtime_error("PALNS route has an invalid request index.");
                }
                const Request& req = data.requests[static_cast<std::size_t>(visit.request_index)];
                if (visit.kind == CpActionKind::Treatment) {
                    if (req.to_id != last_treatment) {
                        treatments.push_back(req.to_id);
                        last_treatment = req.to_id;
                    }
                    treatment_requests.push_back(visit.request_index);
                    continue;
                }
                last_treatment = -1;
                const NodeId next_node = visit.kind == CpActionKind::Pickup ? visit.request_index + 1 : n + visit.request_index + 1;
                if (!seen_nodes.insert(next_node).second) {
                    throw std::runtime_error("PALNS incumbent repeats a service node.");
                }
                flush_to(next_node, req.from_id,
                         visit.kind == CpActionKind::Pickup ? data.time_pickup : data.time_delivery, &visit);
            }
            if (!treatment_requests.empty()) {
                throw std::runtime_error("PALNS route ends with treatment actions.");
            }
            flush_to(graph.end_node_id(), -1, 0.0, nullptr);
            if (!state.is_empty()) {
                throw std::runtime_error("PALNS route does not end empty.");
            }
        }
        if (seen_nodes.size() != 2U * data.requests.size()) {
            throw std::runtime_error("PALNS incumbent does not cover every service node exactly once.");
        }
        const RecoveredSolution recovered = recover_selected_edge_solution(
            data, graph, mapped.active_edge_ids, mapped.total_original_cost, 0.0);
        if (static_cast<int>(recovered.routes.size()) != expected_vehicle_count) {
            throw std::runtime_error("Recovered PALNS route count differs from the vehicle count.");
        }
        for (const RecoveredRouteSolution& r : recovered.routes) {
            if (r.total_time > data.time_limit + 1e-6) {
                throw std::runtime_error("Recovered PALNS route violates the duration limit.");
            }
        }
        mapped.success = true;
    } catch (const std::exception& error) {
        mapped.success = false;
        mapped.error_message = error.what();
        mapped.active_edge_ids.clear();
        mapped.total_duration = 0.0;
        mapped.total_original_cost = 0.0;
    }
    return mapped;
}

}  // namespace spdp

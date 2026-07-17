// recombclust_bigld.hpp — C++ port of the BigLD / CLQD-mod region-detection
// algorithm (Kim et al.; gpart::BigLD / gpart::CLQD), operating on a windows x
// windows correlation matrix produced from ModelsProb.
//
// Boundary (agreed 2026-07-11): everything that BigDataStatMeth can do (HDF5 I/O,
// correlation, colMeans) goes through BDSM in the calling .cpp; the graph
// primitives (coreness, components, max_cliques, graph_from_adjacency_matrix)
// have no BDSM equivalent and are the single documented exception — they are
// delegated to igraph via Rcpp::Environment (same pattern as the legacy
// CGraphUtils.cpp). The pure index/interval logic (thresholding, newSplitCliques,
// constructLDblock, the MWIS in findMaximumIndept, heuristicCLQ) is Eigen/STL on
// the already-loaded, small correlation matrix.
//
// Reference oracle: gpart::CLQD / gpart::BigLD (installed). This port reproduces
// their partitions exactly. NOTE: gpart's LD="r2" branch thresholds abs(cor(geno))
// against CLQcut (NOT cor^2 — the "r2" label is a misnomer); arrange_by_cutoff
// applies abs() so feeding cor(G) with the same CLQcut reproduces its adjacency.

#ifndef RECOMBCLUST_BIGLD_HPP
#define RECOMBCLUST_BIGLD_HPP

#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

namespace recombclust {

// ---------------------------------------------------------------------------
// Small numeric helpers
// ---------------------------------------------------------------------------

// R's type-7 quantile (default of stats::quantile) on an unsorted copy.
inline double quantile7(std::vector<double> x, double p) {
    if (x.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(x.begin(), x.end());
    const std::size_t n = x.size();
    if (n == 1) return x[0];
    double h = (static_cast<double>(n) - 1.0) * p;   // 0-based index position
    std::size_t lo = static_cast<std::size_t>(std::floor(h));
    double frac = h - static_cast<double>(lo);
    if (lo + 1 >= n) return x[n - 1];
    return x[lo] + frac * (x[lo + 1] - x[lo]);
}

// Threshold a correlation matrix into a 0/1 adjacency, mirroring gpart::CLQD:
//   diag -> 0 ; NA -> 0 ; abs(v) < cut -> 0 ; abs(v) >= cut -> 1.
// Vectorised: NaN >= cut is false (IEEE) so missing entries map to 0 for free.
inline Eigen::MatrixXd arrange_by_cutoff(const Eigen::MatrixXd& M, double cut) {
    Eigen::MatrixXd A = (M.array().abs() >= cut).cast<double>();
    A.diagonal().setZero();
    return A;
}

// Symmetric submatrix by index list (0-based) via Eigen 3.4 index slicing.
inline Eigen::MatrixXd submatrix(const Eigen::MatrixXd& M, const std::vector<int>& idx) {
    return M(idx, idx);
}

// ---------------------------------------------------------------------------
// igraph bridge (the documented exception: no BDSM graph primitives)
// ---------------------------------------------------------------------------

inline Rcpp::Environment ig_ns() {
    return Rcpp::Environment::namespace_env("igraph");
}

// Build an undirected igraph from a (weighted, zero-diagonal) adjacency matrix.
inline SEXP ig_graph(const Eigen::MatrixXd& A) {
    Rcpp::Function f = ig_ns()["graph_from_adjacency_matrix"];
    return f(Rcpp::wrap(A), Rcpp::Named("mode", "undirected"),
             Rcpp::Named("weighted", true), Rcpp::Named("diag", false));
}

// k-core number per vertex (unweighted; matches gpart usage).
inline std::vector<int> ig_coreness(SEXP g) {
    Rcpp::Function f = ig_ns()["coreness"];
    Rcpp::NumericVector v = f(g);
    std::vector<int> out(v.size());
    for (int i = 0; i < v.size(); ++i) out[i] = static_cast<int>(std::lround(v[i]));
    return out;
}

// Connected components: membership (0-based here) and component sizes.
inline void ig_components(SEXP g, std::vector<int>& membership, std::vector<int>& csize) {
    Rcpp::Function f = ig_ns()["components"];
    Rcpp::List comp = f(g);
    Rcpp::IntegerVector mem = comp["membership"];   // 1-based component ids
    Rcpp::IntegerVector cs = comp["csize"];
    membership.assign(mem.size(), 0);
    for (int i = 0; i < mem.size(); ++i) membership[i] = mem[i] - 1;
    csize.assign(cs.begin(), cs.end());
}

// Maximal cliques of size >= min_size, as 0-based vertex index lists.
inline std::vector<std::vector<int>> ig_max_cliques(SEXP g, int min_size) {
    Rcpp::Function f = ig_ns()["max_cliques"];
    Rcpp::List cl = f(g, Rcpp::Named("min", min_size));
    std::vector<std::vector<int>> out;
    out.reserve(cl.size());
    for (int i = 0; i < cl.size(); ++i) {
        Rcpp::IntegerVector v = cl[i];
        std::vector<int> c(v.size());
        for (int j = 0; j < v.size(); ++j) c[j] = v[j] - 1;
        out.push_back(std::move(c));
    }
    return out;
}

// All cliques (igraph::cliques) of size >= min_size, 0-based.
inline std::vector<std::vector<int>> ig_cliques(SEXP g, int min_size) {
    Rcpp::Function f = ig_ns()["cliques"];
    Rcpp::List cl = f(g, Rcpp::Named("min", min_size));
    std::vector<std::vector<int>> out;
    out.reserve(cl.size());
    for (int i = 0; i < cl.size(); ++i) {
        Rcpp::IntegerVector v = cl[i];
        std::vector<int> c(v.size());
        for (int j = 0; j < v.size(); ++j) c[j] = v[j] - 1;
        out.push_back(std::move(c));
    }
    return out;
}

// ---------------------------------------------------------------------------
// newSplitCliques — split bp-cliques wherever an internal gap exceeds gapdist.
// Faithful port of gpart:::newSplitCliques. Cliques carry bp values.
// ---------------------------------------------------------------------------
inline std::vector<std::vector<double>>
newSplitCliques(const std::vector<std::vector<double>>& cliques_bp, double gapdist) {
    std::vector<std::vector<double>> nowlist;
    nowlist.reserve(cliques_bp.size());
    for (auto c : cliques_bp) { std::sort(c.begin(), c.end()); nowlist.push_back(std::move(c)); }

    auto max_gap = [](const std::vector<double>& x) -> double {
        if (x.size() < 2) return -std::numeric_limits<double>::infinity();
        double mx = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 1; i < x.size(); ++i) mx = std::max(mx, x[i] - x[i - 1]);
        return mx;
    };

    std::vector<std::vector<double>> fixlist;
    for (;;) {
        std::vector<std::vector<double>> to_split;
        for (auto& x : nowlist) {
            if (max_gap(x) > gapdist) to_split.push_back(x);
            else fixlist.push_back(x);          // need.fix
        }
        if (to_split.empty()) break;

        std::vector<std::vector<double>> newlist;   // rest parts
        std::vector<std::vector<double>> frontparts;
        for (auto& x : to_split) {
            // first gap exceeding gapdist
            std::size_t cut = 0;
            for (std::size_t i = 1; i < x.size(); ++i)
                if (x[i] - x[i - 1] > gapdist) { cut = i; break; }   // cut = # elements in front
            std::vector<double> front(x.begin(), x.begin() + cut);
            std::vector<double> rest(x.begin() + cut, x.end());
            frontparts.push_back(std::move(front));
            newlist.push_back(std::move(rest));
        }
        // front parts with length > 1 are final
        for (auto& f : frontparts) if (f.size() > 1) fixlist.push_back(f);
        // rest parts with length > 1 go around again
        nowlist.clear();
        for (auto& r : newlist) if (r.size() > 1) nowlist.push_back(r);
    }
    return fixlist;
}

// ---------------------------------------------------------------------------
// heuristicCLQ — dense-region extraction (gpart:::heuristicCLQ), used only by
// the near-nonhrst path. Operates on a 0/1 submatrix; returns 0-based indices.
// ---------------------------------------------------------------------------
inline std::vector<int> heuristicCLQ(const Eigen::MatrixXd& subOCM, int hrstParam) {
    const int m = static_cast<int>(subOCM.rows());
    std::vector<double> degs(m);
    for (int i = 0; i < m; ++i) degs[i] = subOCM.row(i).sum();

    double q95 = quantile7(degs, 0.95);
    std::vector<int> top;
    for (int i = 0; i < m; ++i) if (degs[i] >= q95) top.push_back(i);
    std::stable_sort(top.begin(), top.end(),
                     [&](int a, int b) { return degs[a] > degs[b]; });

    auto density = [&](const std::vector<int>& bin) -> double {
        double s = subOCM(bin, bin).sum();
        double L = static_cast<double>(bin.size());
        return s / (L * L - L);
    };

    std::vector<std::vector<int>> basket;
    std::vector<double> dens_all;

    for (int v : top) {
        std::vector<int> v_nbds;
        for (int j = 0; j < m; ++j) if (subOCM(v, j) != 0.0) v_nbds.push_back(j);
        std::vector<int> candibin;
        candibin.push_back(v);
        candibin.insert(candibin.end(), v_nbds.begin(), v_nbds.end());
        double v_density = density(candibin);
        if (v_density >= 0.95) return candibin;

        double new_v_density = 1.0;
        while ((v_density + 1e-4) < new_v_density) {
            // degree of each v_nbds vertex over ALL columns of subOCM
            std::vector<double> newdeg(v_nbds.size());
            double mn = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < v_nbds.size(); ++i) {
                newdeg[i] = subOCM.row(v_nbds[i]).sum();
                mn = std::min(mn, newdeg[i]);
            }
            std::vector<int> kept;
            for (std::size_t i = 0; i < v_nbds.size(); ++i)
                if (newdeg[i] != mn) kept.push_back(v_nbds[i]);
            v_nbds = kept;
            candibin.clear();
            candibin.push_back(v);
            candibin.insert(candibin.end(), v_nbds.begin(), v_nbds.end());
            new_v_density = density(candibin);
            if (new_v_density >= 0.95) return candibin;
            dens_all.push_back(new_v_density);
            basket.push_back(candibin);
            if (static_cast<double>(v_nbds.size()) < (2.0 * hrstParam / 3.0)) break;
        }
    }
    if (basket.empty()) return std::vector<int>();   // defensive
    double mx = *std::max_element(dens_all.begin(), dens_all.end());
    std::vector<int> best_idx;
    for (std::size_t i = 0; i < dens_all.size(); ++i) if (dens_all[i] == mx) best_idx.push_back(i);
    if (best_idx.size() == 1) return basket[best_idx[0]];
    // pick the largest; ties -> first
    int pick = best_idx[0];
    std::size_t bestlen = basket[pick].size();
    for (int bi : best_idx) if (basket[bi].size() > bestlen) { bestlen = basket[bi].size(); pick = bi; }
    return basket[pick];
}

// ---------------------------------------------------------------------------
// CLQD — the core clique-decomposition. Returns a bin vector (0 = unassigned/NA,
// otherwise a 1-based bin id). Faithful port of gpart::CLQD.
//   OCM_raw : correlation matrix (windows x windows), pre-threshold.
//   SNPbps  : bp position per window/SNP (unique).
// ---------------------------------------------------------------------------
inline std::vector<int>
CLQD(const Eigen::MatrixXd& OCM_raw, double CLQcut, double clstgap,
     const std::string& hrstType, int hrstParam, const std::string& CLQmode,
     const std::vector<double>& SNPbps) {

    Eigen::MatrixXd OCM = arrange_by_cutoff(OCM_raw, CLQcut);   // 0/1, diag 0
    const int n = static_cast<int>(OCM.cols());
    std::vector<int> binvector(n, 0);   // 0 == NA
    int binnum = 1;
    if (OCM.sum() == 0.0) return binvector;

    // bp -> global index (positions unique)
    std::map<double, int> bp2idx;
    for (int i = 0; i < n; ++i) bp2idx[SNPbps[i]] = i;

    Eigen::MatrixXd OCM1 = OCM;

    // ---- Heuristic phase ---------------------------------------------------
    if (hrstType == "near-nonhrst") {
        bool heuristic = true;
        while (heuristic) {
            std::vector<int> naidx;
            for (int i = 0; i < n; ++i) if (binvector[i] == 0) naidx.push_back(i);
            // OCM1 is OCM restricted to currently-NA vertices
            OCM1 = submatrix(OCM, naidx);
            SEXP g = ig_graph(OCM1);
            std::vector<int> cores = ig_coreness(g);
            // local_cores: counts of core values >= hrstParam, kept if count >= hrstParam
            std::map<int, int> counts;
            for (int c : cores) if (c >= hrstParam) counts[c]++;
            std::vector<int> local_vals;
            for (auto& kv : counts) if (kv.second >= hrstParam) local_vals.push_back(kv.first);
            if (!local_vals.empty()) {
                int local_hrstParam = *std::max_element(local_vals.begin(), local_vals.end());
                std::vector<int> bothhigh;                 // indices into naidx
                for (std::size_t i = 0; i < cores.size(); ++i)
                    if (cores[i] == local_hrstParam) bothhigh.push_back(static_cast<int>(i));
                Eigen::MatrixXd nowOCM = submatrix(OCM1, bothhigh);
                std::vector<int> hbins = heuristicCLQ(nowOCM, hrstParam);  // idx into bothhigh
                for (int hb : hbins) {
                    int localIdx = bothhigh[hb];      // idx into naidx
                    binvector[naidx[localIdx]] = binnum;
                }
                binnum++;
            } else {
                heuristic = false;
            }
        }
    } else if (hrstType == "fast") {
        bool checkLargest = true;
        // r2Mat and OCM1 shrink together as dense regions are peeled off.
        std::vector<int> rem;                 // remaining global indices (parallels re.SNPbps)
        for (int i = 0; i < n; ++i) rem.push_back(i);
        Eigen::MatrixXd r2Mat = OCM;
        OCM1 = OCM;
        while (checkLargest) {
            SEXP g = ig_graph(r2Mat);
            std::vector<int> membership, csize;
            ig_components(g, membership, csize);
            int componum = 0;                 // 0-based component with max size (first)
            int maxc = -1;
            for (std::size_t i = 0; i < csize.size(); ++i)
                if (csize[i] > maxc) { maxc = csize[i]; componum = static_cast<int>(i); }
            std::vector<int> compov;          // indices into rem
            for (std::size_t i = 0; i < membership.size(); ++i)
                if (membership[i] == componum) compov.push_back(static_cast<int>(i));
            Eigen::MatrixXd compadjM = submatrix(OCM1, compov);
            SEXP cg = ig_graph(compadjM);
            std::vector<int> cgcores = ig_coreness(cg);
            std::vector<double> cgc(cgcores.begin(), cgcores.end());
            double med = quantile7(cgc, 0.5);
            double q75 = quantile7(cgc, 0.75);
            double mx = cgc.empty() ? 0.0 : *std::max_element(cgc.begin(), cgc.end());
            bool fire = ((med > 80 && mx > 100) || (q75 > 100 && mx > 100));
            if (fire) {
                const int R = static_cast<int>(r2Mat.rows());
                std::vector<double> degrees(R);
                for (int i = 0; i < R; ++i) degrees[i] = r2Mat.row(i).sum();
                double thr = std::max(quantile7(degrees, 0.7), 80.0);
                std::vector<int> maxdegvs;
                for (int i = 0; i < R; ++i) if (degrees[i] >= thr) maxdegvs.push_back(i);
                // edge density of each candidate's neighbourhood
                std::vector<double> edgeDens;
                for (int v : maxdegvs) {
                    std::vector<int> nb;
                    for (int j = 0; j < R; ++j) if (r2Mat(v, j) > 0) nb.push_back(j);
                    Eigen::MatrixXd nr = submatrix(r2Mat, nb);
                    double d = static_cast<double>(nr.rows());
                    edgeDens.push_back(nr.sum() / (d * (d - 1)));
                }
                std::vector<int> ord(maxdegvs.size());
                for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = static_cast<int>(i);
                std::stable_sort(ord.begin(), ord.end(),
                                 [&](int a, int b) { return edgeDens[a] > edgeDens[b]; });
                int degv = maxdegvs[ord[0]];
                std::vector<int> nbc;
                for (int j = 0; j < R; ++j) if (r2Mat(degv, j) > 0) nbc.push_back(j);
                // assign these (global via rem) to binnum
                std::vector<int> take_global;
                for (int c : nbc) take_global.push_back(rem[c]);
                for (int gi : take_global) binvector[gi] = binnum;
                binnum++;
                // remove nbc rows/cols from r2Mat, OCM1, rem
                std::vector<bool> drop(R, false);
                for (int c : nbc) drop[c] = true;
                std::vector<int> keep;
                for (int i = 0; i < R; ++i) if (!drop[i]) keep.push_back(i);
                r2Mat = submatrix(r2Mat, keep);
                OCM1 = submatrix(OCM1, keep);
                std::vector<int> newrem;
                for (int i : keep) newrem.push_back(rem[i]);
                rem = newrem;
                checkLargest = (static_cast<int>(rem.size()) >= 500);
            } else {
                checkLargest = false;
            }
        }
    }
    // hrstType == "nonhrst": no heuristic.

    // ---- Split current bins by clstgap ------------------------------------
    bool anyAssigned = false;
    int maxbin = 0;
    for (int b : binvector) { if (b > 0) { anyAssigned = true; maxbin = std::max(maxbin, b); } }
    if (anyAssigned) {
        std::vector<std::vector<double>> binveclist;
        for (int x = 1; x <= maxbin; ++x) {
            std::vector<double> bps;
            for (int i = 0; i < n; ++i) if (binvector[i] == x) bps.push_back(SNPbps[i]);
            binveclist.push_back(std::move(bps));
        }
        std::vector<std::vector<double>> split = newSplitCliques(binveclist, clstgap);
        // keep groups with length > 1, map back to indices, relabel 1..k
        std::fill(binvector.begin(), binvector.end(), 0);
        int lab = 1;
        for (auto& grp : split) {
            if (grp.size() <= 1) continue;
            for (double bp : grp) binvector[bp2idx[bp]] = lab;
            ++lab;
        }
    }

    // ---- Main maximal-clique greedy loop ----------------------------------
    std::vector<int> naidx;
    for (int i = 0; i < n; ++i) if (binvector[i] == 0) naidx.push_back(i);
    Eigen::MatrixXd r2Mat = submatrix(OCM, naidx);
    if (static_cast<int>(naidx.size()) <= 1 || r2Mat.sum() == 0.0) return binvector;

    SEXP g = ig_graph(r2Mat);
    std::vector<std::vector<int>> maxcl = ig_max_cliques(g, 2);   // local idx into naidx
    if (maxcl.empty()) throw std::runtime_error("max.cliques is empty");

    std::vector<double> reSNPbps;                 // bp of currently-NA vertices (fixed)
    for (int i : naidx) reSNPbps.push_back(SNPbps[i]);

    std::vector<std::vector<double>> bp_cliques;
    for (auto& c : maxcl) {
        std::vector<double> v;
        for (int local : c) v.push_back(reSNPbps[local]);
        bp_cliques.push_back(std::move(v));
    }
    std::vector<std::vector<double>> split = newSplitCliques(bp_cliques, clstgap);

    auto sorted_key = [](std::vector<double> v) {
        std::sort(v.begin(), v.end());
        return v;
    };

    for (;;) {
        bool allAssigned = true;
        for (int i = 0; i < n; ++i) if (binvector[i] == 0) { allAssigned = false; break; }
        if (allAssigned) break;
        if (split.empty()) break;

        // density / length per clique
        std::vector<double> dens(split.size());
        for (std::size_t i = 0; i < split.size(); ++i) {
            if (CLQmode == "density") {
                double mn = *std::min_element(split[i].begin(), split[i].end());
                double mxb = *std::max_element(split[i].begin(), split[i].end());
                dens[i] = static_cast<double>(split[i].size()) / ((mxb - mn) / 1000.0);
            } else {
                dens[i] = static_cast<double>(split[i].size());
            }
        }
        double dmax = *std::max_element(dens.begin(), dens.end());
        std::vector<int> cand;
        for (std::size_t i = 0; i < dens.size(); ++i) if (dens[i] == dmax) cand.push_back(static_cast<int>(i));
        int chosen = cand[0];
        std::size_t bestlen = split[chosen].size();
        for (int ci : cand) if (split[ci].size() > bestlen) { bestlen = split[ci].size(); chosen = ci; }
        std::vector<double> max_cluster = split[chosen];

        // remove max_cluster bps from every clique, dedup (stable), keep len > 1
        std::vector<double> mc = sorted_key(max_cluster);
        std::vector<std::vector<double>> next;
        std::vector<std::vector<double>> seen;
        for (auto& c : split) {
            std::vector<double> diff;
            for (double bp : c)
                if (!std::binary_search(mc.begin(), mc.end(), bp)) diff.push_back(bp);
            if (diff.size() <= 1) continue;
            std::vector<double> key = sorted_key(diff);
            bool dup = false;
            for (auto& s : seen) if (s == key) { dup = true; break; }
            if (dup) continue;
            seen.push_back(key);
            next.push_back(diff);
        }

        for (double bp : max_cluster) binvector[bp2idx[bp]] = binnum;
        ++binnum;
        split = next;

        if (static_cast<int>(reSNPbps.size()) < 2) break;
        if (split.empty()) break;
    }
    return binvector;
}

// ---------------------------------------------------------------------------
// intervalCliqueList — gpart:::intervalCliqueList. clstlist holds integer index
// sets (their VALUES matter; they need not be 0/1-based). Returns rows
// {min, max, weight} of candidate intervals (deduped by range, max weight).
// ---------------------------------------------------------------------------
inline std::vector<std::array<int, 3>>
intervalCliqueList(std::vector<std::vector<int>> clstlist) {
    auto vmin = [](const std::vector<int>& v) { return *std::min_element(v.begin(), v.end()); };
    auto vmax = [](const std::vector<int>& v) { return *std::max_element(v.begin(), v.end()); };

    std::stable_sort(clstlist.begin(), clstlist.end(),
                     [&](const std::vector<int>& a, const std::vector<int>& b) {
                         return vmin(a) < vmin(b);
                     });
    const int IM = static_cast<int>(clstlist.size());

    Eigen::MatrixXd adj = Eigen::MatrixXd::Zero(IM, IM);
    for (int i = 0; i < IM - 1; ++i)
        for (int j = i + 1; j < IM; ++j)
            if (vmax(clstlist[i]) > vmin(clstlist[j])) { adj(i, j) = 1.0; adj(j, i) = 1.0; }

    SEXP g = ig_graph(adj);
    std::vector<int> cores = ig_coreness(g);
    int maxcore = cores.empty() ? 0 : *std::max_element(cores.begin(), cores.end());
    std::vector<std::vector<int>> cliques =
        (maxcore > 10) ? ig_max_cliques(g, 1) : ig_cliques(g, 1);

    // order cliques by their minimum vertex index
    std::stable_sort(cliques.begin(), cliques.end(),
                     [](const std::vector<int>& a, const std::vector<int>& b) {
                         return *std::min_element(a.begin(), a.end()) <
                                *std::min_element(b.begin(), b.end());
                     });

    // For each clique: union member clusters' values -> sorted unique interval.
    std::vector<int> mins, maxs, weights;
    for (auto& cl : cliques) {
        std::vector<int> vals;
        for (int v : cl) vals.insert(vals.end(), clstlist[v].begin(), clstlist[v].end());
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
        mins.push_back(vals.front());
        maxs.push_back(vals.back());
        weights.push_back(static_cast<int>(vals.size()));
    }

    // Dedup by (min,max) in first-occurrence order, keeping max weight.
    std::vector<std::array<int, 3>> res;
    std::vector<std::pair<int, int>> seen;
    for (std::size_t i = 0; i < mins.size(); ++i) {
        std::pair<int, int> key(mins[i], maxs[i]);
        if (std::find(seen.begin(), seen.end(), key) != seen.end()) continue;
        seen.push_back(key);
        int w = -1;
        for (std::size_t k = 0; k < mins.size(); ++k)
            if (mins[k] == key.first && maxs[k] == key.second) w = std::max(w, weights[k]);
        res.push_back({key.first, key.second, w});
    }
    return res;
}

// ---------------------------------------------------------------------------
// findMaximumIndept — gpart:::findMaximumIndept. Maximum-weight independent set
// over intervals via a DAG longest-path, replicated in gpart's index order.
// ranges: {start,end} rows; weights: per-interval weight. Returns the 0-based
// interval indices forming the chosen set.
// ---------------------------------------------------------------------------
inline std::vector<int>
findMaximumIndept(const std::vector<std::array<int, 2>>& ranges,
                  const std::vector<int>& weights) {
    const int n = static_cast<int>(weights.size());
    std::vector<std::vector<int>> pre(n);           // predecessors (empty == NA)
    for (int i = 0; i < n; ++i) {
        int nowstart = ranges[i][0];
        for (int k = 0; k < n; ++k)
            if (ranges[k][1] < nowstart) pre[i].push_back(k);
    }
    std::vector<int> sources, not_s;
    for (int i = 0; i < n; ++i) (pre[i].empty() ? sources : not_s).push_back(i);

    std::vector<int> indept;
    if (static_cast<int>(sources.size()) < n) {
        // keep only immediate predecessors (drop predecessors-of-predecessors)
        for (int i : not_s) {
            std::vector<int> prepre;
            for (int p : pre[i]) prepre.insert(prepre.end(), pre[p].begin(), pre[p].end());
            std::sort(prepre.begin(), prepre.end());
            prepre.erase(std::unique(prepre.begin(), prepre.end()), prepre.end());
            std::vector<int> keep;
            for (int p : pre[i])
                if (!std::binary_search(prepre.begin(), prepre.end(), p)) keep.push_back(p);
            pre[i] = keep;
        }
        std::vector<double> route(n, 0.0);
        std::vector<int> pointer(n, -1);            // -1 == NA
        for (int s : sources) { route[s] = weights[s]; pointer[s] = -1; }

        for (int i : not_s) {                       // ascending index (gpart order)
            const std::vector<int>& mp = pre[i];
            int best = mp[0];
            double bestw = route[mp[0]];
            for (int p : mp) if (route[p] > bestw) { bestw = route[p]; best = p; }
            route[i] = weights[i] + bestw;
            pointer[i] = best;
        }
        int start_itv = 0;
        double mx = route[0];
        for (int i = 1; i < n; ++i) if (route[i] > mx) { mx = route[i]; start_itv = i; }
        // walk the pointer chain back
        std::vector<int> chain;
        int predecessor = pointer[start_itv];
        chain.push_back(predecessor);
        chain.push_back(start_itv);
        while (predecessor != -1) {
            predecessor = pointer[predecessor];
            chain.insert(chain.begin(), predecessor);
        }
        for (int c : chain) if (c != -1) indept.push_back(c);
    } else {
        int mx = *std::max_element(weights.begin(), weights.end());
        for (int i = 0; i < n; ++i) if (weights[i] == mx) indept.push_back(i);
    }
    return indept;
}

// ---------------------------------------------------------------------------
// constructLDblock — gpart:::constructLDblock. clstlist holds index sets; returns
// LD blocks as {start,end} index ranges (unordered, as gpart).
// ---------------------------------------------------------------------------
inline std::vector<std::array<int, 2>>
constructLDblock(std::vector<std::vector<int>> clstlist) {
    auto vmin = [](const std::vector<int>& v) { return *std::min_element(v.begin(), v.end()); };
    auto vmax = [](const std::vector<int>& v) { return *std::max_element(v.begin(), v.end()); };

    std::vector<std::array<int, 2>> total;
    while (!clstlist.empty()) {
        if (clstlist.size() == 1) {
            total.push_back({vmin(clstlist[0]), vmax(clstlist[0])});
            break;
        }
        std::vector<std::array<int, 3>> candi = intervalCliqueList(clstlist);
        std::vector<std::array<int, 2>> intervals;
        std::vector<int> weight;
        for (auto& c : candi) { intervals.push_back({c[0], c[1]}); weight.push_back(c[2]); }
        std::vector<int> mwis = findMaximumIndept(intervals, weight);

        std::vector<std::array<int, 2>> subblocks;
        for (int i : mwis) subblocks.push_back(intervals[i]);
        for (auto& b : subblocks) total.push_back(b);

        // takenSNPs = every index in [min,max] of each chosen block
        std::vector<int> takenv;
        for (auto& b : subblocks) for (int x = b[0]; x <= b[1]; ++x) takenv.push_back(x);
        std::sort(takenv.begin(), takenv.end());
        takenv.erase(std::unique(takenv.begin(), takenv.end()), takenv.end());
        auto taken = [&](int x) { return std::binary_search(takenv.begin(), takenv.end(), x); };

        // remove taken from every cluster, drop empties
        std::vector<std::vector<int>> nxt;
        for (auto& c : clstlist) {
            std::vector<int> d;
            for (int e : c) if (!taken(e)) d.push_back(e);
            if (!d.empty()) nxt.push_back(std::move(d));
        }
        clstlist = std::move(nxt);
        if (clstlist.empty()) break;

        // re-fragment clusters straddling a taken region until stable
        for (;;) {
            std::vector<std::vector<int>> adding;
            for (auto& nowbin : clstlist) {
                int lo = vmin(nowbin), hi = vmax(nowbin);
                int minInter = 0, maxInter = 0; bool hasInter = false;
                for (int x = lo; x <= hi; ++x)
                    if (taken(x)) {
                        if (!hasInter) { minInter = x; hasInter = true; }
                        maxInter = x;
                    }
                if (hasInter) {
                    std::vector<int> part1, part2;
                    for (int e : nowbin) {
                        if (e < minInter) part1.push_back(e);
                        else if (e > maxInter) part2.push_back(e);
                    }
                    nowbin = part1;
                    adding.push_back(std::move(part2));
                }
            }
            for (auto& a : adding) clstlist.push_back(std::move(a));
            std::vector<std::vector<int>> keep;
            for (auto& c : clstlist) if (c.size() > 1) keep.push_back(std::move(c));
            clstlist = std::move(keep);
            std::stable_sort(clstlist.begin(), clstlist.end(),
                             [&](const std::vector<int>& a, const std::vector<int>& b) {
                                 return vmin(a) < vmin(b);
                             });
            if (adding.empty()) break;
        }
    }
    return total;
}

// ---------------------------------------------------------------------------
// BigLD — CLQD -> clstlist -> constructLDblock -> LD blocks, mirroring the
// recombClust BigLD wrapper. win_start / win_end are the bp span of each window
// (SNPbps for CLQD = win_start). Returns blocks as {start_idx, end_idx} (0-based)
// sorted by start; the caller maps indices to bp via win_start/win_end.
// ---------------------------------------------------------------------------
inline std::vector<std::array<int, 2>>
BigLD(const Eigen::MatrixXd& OCM, double CLQcut, double clstgap,
      const std::string& hrstType, int hrstParam, const std::string& CLQmode,
      const std::vector<double>& win_start) {
    std::vector<int> binvector =
        CLQD(OCM, CLQcut, clstgap, hrstType, hrstParam, CLQmode, win_start);
    int maxbin = 0;
    for (int b : binvector) maxbin = std::max(maxbin, b);
    if (maxbin == 0) return std::vector<std::array<int, 2>>();

    // clstlist: 1-based index sets per bin (to match gpart's which()), sorted.
    std::vector<std::vector<int>> clstlist;
    for (int x = 1; x <= maxbin; ++x) {
        std::vector<int> idx;
        for (int i = 0; i < static_cast<int>(binvector.size()); ++i)
            if (binvector[i] == x) idx.push_back(i + 1);
        if (!idx.empty()) { std::sort(idx.begin(), idx.end()); clstlist.push_back(std::move(idx)); }
    }
    std::stable_sort(clstlist.begin(), clstlist.end(),
                     [](const std::vector<int>& a, const std::vector<int>& b) {
                         return a.front() < b.front();
                     });

    std::vector<std::array<int, 2>> blocks = constructLDblock(clstlist);
    // back to 0-based indices, sort by start
    for (auto& b : blocks) { b[0] -= 1; b[1] -= 1; }
    std::stable_sort(blocks.begin(), blocks.end(),
                     [](const std::array<int, 2>& a, const std::array<int, 2>& b) {
                         return a[0] < b[0];
                     });
    return blocks;
}

}  // namespace recombclust

#endif  // RECOMBCLUST_BIGLD_HPP

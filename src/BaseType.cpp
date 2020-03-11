#include "BaseType.h"
#define FMT_HEADER_ONLY
#include "fmt/format.h"

BaseType::BaseType(BaseV bases_, BaseV quals_, int8_t ref, double minaf) : bases(bases_), quals(quals_), ref_base(ref), min_af(minaf), nind(bases.size()), ind_allele_likelihood(nind * NTYPE), init_allele_freq(NTYPE)
{
    var_qual = 0;
    depth_total = 0;
    depth = { {0, 0},{1, 0},{2, 0},{3, 0} };
    for (int32_t i = 0; i < nind; ++i) {
        for (int j = 0; j < NTYPE; ++j) {
            if (bases[i] == BASE[j]) {
                ind_allele_likelihood[i * NTYPE + j] = 1.0 - exp(MLN10TO10 * quals[i]);
            } else {
                ind_allele_likelihood[i * NTYPE + j] = exp(MLN10TO10 * quals[i]) / 3.0;
            }
        }
        depth[bases[i]] += 1;
    }
    for (DepM::iterator it = depth.begin(); it != depth.end(); ++it) {
        depth_total += it->second;
    }
}

void BaseType::SetAlleleFreq(const BaseV& bases)
{
    int32_t depth_sum = 0;
    for (auto b : bases) {
        depth_sum += depth[b];
    }
    for (int j = 0; j < NTYPE; ++j) {
        init_allele_freq[j] = 0;
    }
    if (depth_sum > 0) {
        for (auto b : bases) {
            init_allele_freq[b] = static_cast<double>(depth[b]) / depth_sum;
        }
    }
}

void BaseType::UpdateF(const BaseV& bases, CombV& bc, ProbV& lr, FreqV& bp, int32_t k)
{
    ProbV marginal_likelihood(nind), expect_allele_prob(NTYPE);
    double freq_sum, likelihood_sum;
    double epsilon = 0.001;
    int iter_num = 100;
    ProbV expect_prob;
    bc.clear(); lr.clear(); bp.clear();
    combs_(bases, bc, k);
    for (auto const& b: bc) {
        SetAlleleFreq(b);
        freq_sum = 0;
        for (int i = 0; i < NTYPE; ++i) freq_sum += init_allele_freq[i];
        if (freq_sum == 0) continue;  // skip coverage = 0, this may be redundant but it's ok;
        // run EM
        EM(init_allele_freq, ind_allele_likelihood, marginal_likelihood, expect_allele_prob, nind, NTYPE, iter_num, epsilon);
        likelihood_sum = 0;
        for (int32_t i = 0; i < nind; ++i) {
            likelihood_sum += std::log(marginal_likelihood[i]);
            // reset array elements to 0
            marginal_likelihood[i] = 0;
        }
        expect_prob.clear();
        for (int i = 0; i < NTYPE; ++i) {
            expect_prob.push_back(expect_allele_prob[i]);
            expect_allele_prob[i] = 0;
        }
        lr.push_back(likelihood_sum);
        bp.push_back(expect_prob);
    }
}

bool BaseType::LRT()
{
    if (depth_total == 0) return false;
    bases.clear();
    for (auto b : base_comb) {
        // filter bases by count freqence >= min_af
        if ((depth[b]/depth_total) >= min_af) {
            bases.push_back(b);
        }
    }
    int32_t n = bases.size();
    if (n == 0) return false;
    CombV bc;
    FreqV bp;
    ProbV lr_null, lrt_chi;
    UpdateF(bases, bc, lr_null, bp, n);
    ProbV base_frq = bp[0];
    double lr_alt_t = lr_null[0];
    double chi_sqrt_t = 0.0;
    size_t i_min;
    for (int32_t k = n - 1; k > 0; --k) {
        UpdateF(bases, bc, lr_null, bp, k);
        lrt_chi.clear();
        for (auto & lr_null_t: lr_null) {
            lrt_chi.push_back(2.0 * (lr_alt_t - lr_null_t));
        }
        i_min = std::min_element(lrt_chi.begin(), lrt_chi.end()) - lrt_chi.begin();
        lr_alt_t = lr_null[i_min];
        chi_sqrt_t = lrt_chi[i_min];
        if (chi_sqrt_t < LRT_THRESHOLD) {
            // Take the null hypothesis and continue
            bases = bc[i_min];
            base_frq = bp[i_min];
        } else {
            // Take the alternate hypothesis
            break;
        }
    }
    for (auto b: bases) {
        if (b != ref_base) {
            alt_bases.push_back(b);
            af_lrt.insert({b, base_frq[b]});
        }
    }
    double r, chi_prob;
    if (alt_bases.size() > 0) {
        r = depth[bases[0]] / depth_total;
        if (bases.size() == 1 && depth_total > 10 && r > 0.5) {
            var_qual = 5000.0;  // mono-allelelic
        } else {
            if (chi_sqrt_t <= 0) {
                var_qual = 0.0;
                return true;
            }
            chi_prob = chisf(chi_sqrt_t, 1.0);  // may be nan value;
            if (chi_prob) {
                var_qual = -10 * log10(chi_prob);
            } else {
                var_qual = 10000;
            }
            if (var_qual == 0) var_qual = 0.0;  // output -0.0 to 0;
        }
        return true;
    } else {
        return false;
    }
}

String WriteVcf(const BaseType& bt, const String& chr, int32_t pos, int8_t ref_base, const AlleleInfoVector& aiv, const DepM& idx, InfoM& info, int32_t N)
{
    robin_hood::unordered_map<uint8_t, String> alt_gt;
    String gt, samgt;
    for (size_t i = 0; i < bt.alt_bases.size(); ++i) {
        gt = fmt::format("./{}", i+1);
        alt_gt.insert({bt.alt_bases[i], gt});
    }
    Stat st;
    ProbV ref_quals, ref_mapqs, ref_rprs;
    ProbV alt_quals, alt_mapqs, alt_rprs;
    for (int32_t i = 0; i < N; ++i) {
        if (idx.count(i) == 0) {
            samgt += "./.\t";
        } else {
            auto const& a = aiv[idx.at(i)];
            if (alt_gt.count(a.base) == 0) alt_gt.insert({a.base, "./."});
            if (a.base == ref_base) {
                gt = "0/.";
            } else {
                gt = alt_gt[a.base];
            }
            samgt += fmt::format("{}:{}:{}:{:.6f}\t", gt, BASE2CHAR[a.base], STRAND[a.strand], 1 - exp(MLN10TO10 * a.qual));
            if (a.is_indel == 1 || a.base == 4) continue;
            if (a.base == ref_base) {
                ref_quals.push_back(a.qual);
                ref_mapqs.push_back(a.mapq);
                ref_rprs.push_back(a.rpr);
            } else if (std::find(bt.alt_bases.begin(), bt.alt_bases.end(), a.base) != bt.alt_bases.end()) {
                alt_quals.push_back(a.qual);
                alt_mapqs.push_back(a.mapq);
                alt_rprs.push_back(a.rpr);
            }
            if (a.strand == 1) {
                if (a.base == ref_base) {
                    st.ref_fwd += 1;
                } else if (std::find(bt.alt_bases.begin(), bt.alt_bases.end(), a.base) != bt.alt_bases.end()) {
                    st.alt_fwd += 1;
                }
            } else if (a.strand == 0) {
                if (a.base == ref_base) {
                    st.ref_rev += 1;
                } else if (std::find(bt.alt_bases.begin(), bt.alt_bases.end(), a.base) != bt.alt_bases.end()) {
                    st.alt_rev += 1;
                }
            }
        }
    }
    st.phred_mapq = RankSumTest(ref_mapqs, alt_mapqs);
    st.phred_qual = RankSumTest(ref_quals, alt_quals);
    st.phred_rpr = RankSumTest(ref_rprs, alt_rprs);
    st.fs = bt_fisher_exact(st.ref_fwd, st.ref_rev, st.alt_fwd, st.alt_rev);
    if (st.alt_fwd * st.ref_rev > 0) {
        st.sor = static_cast<double>(st.ref_fwd * st.alt_rev) / (st.ref_rev * st.alt_fwd);
    } else {
        st.sor = 10000.0;
    }
    double ad_sum = 0;
    String ac, af, caf, alt;
    for (auto b : bt.alt_bases) {
        ad_sum += bt.depth.at(b);
        alt += fmt::format("{},", BASE2CHAR[b]);
        ac += fmt::format("{},", bt.depth.at(b));
        af += fmt::format("{:.6f},", bt.af_lrt.at(b));
        caf += fmt::format("{:.6f},", bt.depth.at(b) / bt.depth_total);
    }
    alt.pop_back(); samgt.pop_back();
    ac.pop_back(); info.insert({"CM_AC", ac});
    af.pop_back(); info.insert({"CM_AF", af});
    caf.pop_back(); info.insert({"CM_CAF", caf});
    info.insert({"QD", fmt::format("{:.3f}", bt.var_qual/ad_sum)});
    info.insert({"CM_DP", fmt::format("{:.0f}", bt.depth_total)});
    info.insert({"MQRankSum", fmt::format("{:.3f}", st.phred_mapq)});
    info.insert({"ReadPosRankSum", fmt::format("{:.3f}", st.phred_rpr)});
    info.insert({"BaseQRankSum", fmt::format("{:.3f}", st.phred_qual)});
    info.insert({"FS", fmt::format("{:.3f}", st.fs)});
    info.insert({"SOR", fmt::format("{:.3f}", st.sor)});
    info.insert({"SB_REF", fmt::format("{},{}", st.ref_fwd, st.ref_rev)});
    info.insert({"SB_ALT", fmt::format("{},{}", st.alt_fwd, st.alt_rev)});
    String qt;
    if (bt.var_qual > QUAL_THRESHOLD) {
        qt = ".";
    } else {
        qt = "LowQual";
    }
    String out = fmt::format("{}\t{}\t.\t{}\t{}\t{:.2f}\t{}\t", chr, pos, BASE2CHAR[ref_base], alt, bt.var_qual, qt);
    for (InfoM::iterator it = info.begin(); it != info.end(); ++it) {
        out += it->first + "=" + it->second + ";";
    }
    out.pop_back();
    out += "\tGT:AB:SO:BP\t" + samgt + "\n";

    return out;
}


void combs_(const BaseV& bases, CombV& comb_v, int32_t k)
{
    // k <= n
    int32_t n = bases.size();
    std::string bitmask(k, 1); // K leading 1's
    bitmask.resize(n, 0); // N-K trailing 0's

    BaseV bv;
    comb_v.clear();
    bv.reserve(k);
    do {
        for (int32_t i = 0; i < n; ++i) // [0..N-1] integers
        {
            if (bitmask[i]) bv.push_back(bases[i]);
        }
        comb_v.push_back(bv);
        bv.clear();
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

#include "BaseType.h"

BaseType::BaseType(BaseV base, BaseV qual, int ref, double minaf) : bases(base), quals(qual), ref_base(ref), min_af(minaf)
{
    nind = bases.size();
    var_qual = 0;
    ind_allele_likelihood = new double[nind * NTYPE]; // @WATCHOUT
    for (int32_t i = 0; i < nind; ++i) {
        for (int j = 0; j < NTYPE; ++j) {
            if (bases[i] == BASE[j]) {
                ind_allele_likelihood[i * NTYPE + j] = 1.0 - exp(MLN10TO10 * quals[i]);
            } else {
                ind_allele_likelihood[i * NTYPE + j] = exp(MLN10TO10 * quals[i]) / 3.0;
            }
        }
        if (depth.count(bases[i]) == 1) depth[bases[i]] += 1;
    }
    for (DepM::iterator it = depth.begin(); it != depth.end(); ++it) {
        depth_total += it->second;
    }
}

void BaseType::SetAlleleFreq(const BaseV& bases) {
    double depth_sum = 0;
    for (auto const& b : bases) {
        depth_sum += depth[b];
    }
    for (int j = 0; j < NTYPE; ++j) {
        init_allele_freq[j] = 0;
    }
    if (depth_sum > 0) {
        for (auto const& b : bases) {
            init_allele_freq[b] = depth[b] / depth_sum;
        }
    }
}

void BaseType::UpdateF(const BaseV& bases, CombV& bc, ProbV& lr, FreqV& bp, int32_t k) {
    // REMINDME: consider comb_v.reserve
    combs(bases, bc, k);
    double *marginal_likelihood = new double[nind](); // @WATCHOUT
    double *expect_allele_prob = new double[NTYPE]();
    double freq_sum = 0;
    double likelihood_sum = 0;
    double epsilon = 0.001;
    int iter_num = 100;
    ProbV expect_prob;
    lr.clear();
    bp.clear();
    for (auto const& b: bc) {
        SetAlleleFreq(b);
        freq_sum = 0;
        for (int i = 0; i < NTYPE; ++i) freq_sum += init_allele_freq[i];
        if (freq_sum == 0) continue;
        // run EM
        EM(init_allele_freq, ind_allele_likelihood, marginal_likelihood, expect_allele_prob, nind, NTYPE, iter_num, epsilon);
        likelihood_sum = 0;
        for (int32_t i = 0; i < nind; ++i) {
            likelihood_sum += log(marginal_likelihood[i]);
            // reset array elements to 0
            marginal_likelihood[i] = 0;
        }
        expect_prob.assign(expect_allele_prob, expect_allele_prob + NTYPE);
        lr.push_back(likelihood_sum);
        bp.push_back(expect_prob);
        std::fill(expect_allele_prob, expect_allele_prob + NTYPE, 0);
    }

    delete []marginal_likelihood;
    delete []expect_allele_prob;
}

void BaseType::LRT(){
    if (depth_total == 0) return;
    if (nind == 0) return;
    CombV bc;
    FreqV bp;
    ProbV lr_null;
    ProbV lrt_chi;
    ProbV base_frq;
    UpdateF(bases, bc, lr_null, bp, nind);
    base_frq = bp[0];
    auto lr_alt_t = lr_null[0];
    double chi_sqrt_t = 0;
    int i_min;
    for (int32_t k = bases.size() - 1; k > 0; --k) {
        UpdateF(bases, bc, lr_null, bp, k);
        lrt_chi.clear();
        for (auto& lr_null_t: lr_null) {
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
    for (auto& b: bases) {
        if (b != ref_base) {
            alt_bases.push_back(b);
            af_lrt.push_back(base_frq[b]);
            // af_lrt.insert({b, base_frq[b]});
        }
    }
    double r, chi_prob;
    if (alt_bases.size() > 0) {
        r = depth[bases[0]] / depth_total;
        if (bases.size() == 1 && depth_total > 10 && r > 0.5) {
            var_qual = 5000.0;  // mono-allelelic
        } else {
            chi_prob = chisf(chi_sqrt_t, 1.0);
            if (chi_prob) {
                var_qual = -10 * log10(chi_prob);
            } else {
                var_qual = 10000.0;
            }
        }
    }
}
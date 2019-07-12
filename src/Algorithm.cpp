#include "Algorithm.h"

double chisf(double x, double k)
{
    // sf = 1 - cdf
    return kf_gammaq(k/2.0, x/2.0);
}

double normsf(double x)
{
    // mean = 0, sd = 1
    return 1 - (1 + erf(x/std::sqrt(2.0)))/2;
}

double bt_fisher_exact(int n11, int n12, int n21, int n22)
{
    double left_p, right_p, twoside_p;
    kt_fisher_exact(n11, n12, n21, n22, &left_p, &right_p, &twoside_p);
    double p = -10 * std::log10(twoside_p);
    if (std::isinf(p)) p = 10000.0;
    else if (p == 0) p = 0.0;

    return p;
}

static double rankR1(const std::vector<double>& x, size_t n1)
{
    // merge, sort and keep track of index
    std::vector<size_t> idx = BaseVar::sortidx(x);
    // average method
    int32_t i, j, k = 0, n = 0, s = idx.size();
    size_t id1, id2;
    double avg, r1 = 0.0;
    for (i = 0; i < s; ++i) {
        id1 = idx[i]; id2 = idx[i+1];
        if (i + 1 < s && x[id1] == x[id2]) {
            k += i + 1;
            n++;
        } else {
            if (k > 0) {
                k += i + 1;
                avg = static_cast<double>(k) / (n + 1);
                for (j = i; i-n-j <= 0; j--) if (idx[j] < n1) r1 += avg;
                k = 0; n = 0;
            } else {
                if (id1 < n1) r1 += i + 1;
            }
        }
    }

    return r1;
}

double RankSumTest(std::vector<double>& x, std::vector<double>& y)
{
    size_t n1 = x.size(), n2 = y.size();
    x.insert(x.end(), y.begin(), y.end());
    double r1 = rankR1(x, n1);
    double expected = (double)(n1 * (n1 + n2 + 1)) / 2.0;
    double z = (r1 - expected) / std::sqrt(static_cast<double>(n1*n2*(n1+n2+1))/12.0);
    double p = -10 * std::log10(2 * normsf(std::abs(z))); // phred score
    if (std::isinf(p)) p = 10000.0;
    else if (!(p != 0)) p = 0.0;

    return p;
}

static void singleEM(const std::vector<double>& allele_freq, const std::vector<double>& ind_allele_likelihood, std::vector<double>& marginal_likelihood, std::vector<double>& expect_allele_prob, int32_t nsample, int ntype)
{
    std::vector<double> likelihood(ntype), ind_allele_prob(ntype * nsample);
    // step E
    int i, j;
    for (i = 0; i < nsample; ++i) {
        for (j = 0; j < ntype; ++j) {
            likelihood[j] = allele_freq[j] * ind_allele_likelihood[i * ntype + j];
            marginal_likelihood[i] += likelihood[j];
        }
        for (j = 0; j < ntype; ++j) {
            // col major may be fast for step M
            // need to deal with marginal_likelihood[i] is close to zero
            ind_allele_prob[j * nsample + i] = likelihood[j] / marginal_likelihood[i];
        }
    }
    // step M
    for (j = 0; j < ntype; ++j){
        for (i = 0; i < nsample; ++i){
            expect_allele_prob[j] += ind_allele_prob[j * nsample + i];
        }
        expect_allele_prob[j] = expect_allele_prob[j] / nsample;
    }
    return;
}

static inline void update_allele_freq(std::vector<double>& allele_freq, std::vector<double>& expect_allele_prob, int ntype)
{
    for(int j = 0; j < ntype; ++j) {
        allele_freq[j] = expect_allele_prob[j];
        expect_allele_prob[j] = 0.0;
    }
}

static double delta_bylog(std::vector<double>& bf, std::vector<double>& af, int nsample)
{
    double delta = 0.0;
    for(int i = 0; i < nsample; ++i) {
        // need to deal with log(0) == inf;
        delta += std::abs(std::log(af[i]) - std::log(bf[i]));
        bf[i] = af[i];
        af[i] = 0.0;
    }
    return delta;
}

void EM(std::vector<double>& init_allele_freq, const std::vector<double>& ind_allele_likelihood, std::vector<double>& marginal_likelihood, std::vector<double>& expect_allele_prob, int32_t nsample, int ntype, int iter_num, double epsilon)
{
    std::vector<double> af_marginal_likelihood(nsample);
    singleEM(init_allele_freq, ind_allele_likelihood, marginal_likelihood, expect_allele_prob, nsample, ntype);
    double delta;
    for(int i = 0; i < iter_num; ++i){
        update_allele_freq(init_allele_freq, expect_allele_prob, ntype);
        singleEM(init_allele_freq, ind_allele_likelihood, af_marginal_likelihood, expect_allele_prob, nsample, ntype);
        delta = delta_bylog(marginal_likelihood, af_marginal_likelihood, nsample);
        if(delta < epsilon){
            break;
        }
    }

    return;
}

